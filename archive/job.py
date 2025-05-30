import os
import sys  # noqa: F401
from math import floor, ceil, sqrt, log  
import numpy as np
from ase.io import read, write  # noqa: F401
#from ase.visualize import view
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.optimize import QuasiNewton
from pathlib import Path

# Set directories for importable scripts and templates 
cwd = Path.cwd()
project_dir = cwd.parent.parent.parent.parent
sys.path.insert(1, project_dir.as_posix())

# Import custom scripts
from lib.pgroup.requeue import ReQueue  # noqa: E402

def write_energy(format_string, *args):
    """ Write the energy at a certain step in the calculation.
    'status' should indicate, in some way, "Initial" or "Final".
     Example usage:
    write_energy('{} Energy: {}, {} Forces: {}', 'Final', 'Final', 42, 3.14)
    """
    # Create a Path object for the specified filename
    file_path = Path('energy.log')

    # Format passed floats
    formatted_args = []
    for arg in args:
        if isinstance(arg, float):
            formatted_args.append(f'{arg:.4f}')
        else:
            formatted_args.append(str(arg))

    # Format the string with the provided arguments
    formatted_line = format_string.format(*formatted_args) + '\n'
        
    # Append the formatted line to the file
    with file_path.open(mode='a', encoding='utf-8') as file:
        file.write(formatted_line)


# Check if vasp path set correctly, if not, exit early
vasp_path = os.environ.get('EBROOTVASP') #EBROOTVASP, Niagara = SCINET_VASP_ROOT
if not isinstance(vasp_path, str):  
    exit()

num_cores = int(os.environ.get('SLURM_NPROCS'))
npar_setting = 2**ceil(log(floor(sqrt(num_cores)), 2)) # NPAR must be a divisor of the number of cores. Assumed to be a power of 2.

# Read system from prior opt, clear old calc and constraints
num_OH = 3
system = read(f'Pt_{num_OH}_OH_global_min.traj')
system.calc = None
system.constraints = None

# Fix bottom two layers, check structure
constrain_layers = [atom.tag > 2 for atom in system]
system.set_constraint(FixAtoms(mask=constrain_layers))
#view(system)

# Construct the calculator with necessary parameters
calc = Vasp(prec = 'Accurate',
            xc = 'PBE',
            encut = 650, #Planewave cutoff to be set in loop
            ispin = 2, # 1 = no spin-polarization, default
            ibrion = -1, # -1 = no Vasp relaxation. Defult = -1 iff nsw = 0
            nsw = 0, # Max no. relaxation steps. This is a default
            kpts = [6, 6, 1],
            ismear = 0, #https://www.vasp.at/wiki/index.php/ISMEAR
            sigma = 0.05,
            #potim = 0.5, #time step/step width for md and relaxation
            isif = 2, #Calc stress tensors
            npar = npar_setting # Band parallelization recommend to set NPAR on these machines to √# of cores
)

#Set the trajectory, log file, and calc
calc.set(txt=f'Pt_{num_OH}_OH_final_opt.txt')
system.calc = calc
opt_traj = f'Pt_{num_OH}_OH_final_opt.traj'
opt_log = f'Pt_{num_OH}_OH_final_opt.log'
opt_hessian='qn.json'
opt = QuasiNewton(system,
                  trajectory=opt_traj,
                  logfile=opt_log,
                  restart=opt_hessian)
system.calc = calc
    
#Get and Store Initial Energy
e = system.get_potential_energy()
print(f'Num OH: {num_OH}. Energy: {e:.2f} eV.')
write_energy('Num OH: {}, Initial Energy: {} eV', num_OH, e)

# Run relaxation with QuasiNewton and Requeuing
requeue = ReQueue(maxtime=23., checktime=0.5, log='requeue.log')
# Start Compute Intensive Portion of Job
status = requeue(opt.run, fmax=0.01)

# Get status, requeue if timed out, get final Energy if done
if status == 'time_elapsed':
	os.system("sbatch --begin=now+30minutes vaspJob.sh")
elif status == 'job_completed':
	#Get and Store Final Energy, Max Force
	e = system.get_potential_energy()
	f = system.get_forces()
	f_max = np.max(np.linalg.norm(f, axis=1))
	write_energy('Num OH: {}, Final Energy: {} eV, Final Max Force: {} eV/A', num_OH, e, f_max)
	print(f'Num OH: {num_OH}, Energy: {e:.2f} eV. Max Force: {f_max:.4f}')