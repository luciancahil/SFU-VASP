#! /usr/bin/env python

from ase.io import read, Trajectory, write
from ase.build import fcc111, root_surface,add_adsorbate, add_vacuum
from ase import Atoms
import ase.calculators.vasp as vasp_calculator
from ase.optimize import BFGS
from ase.build import surface
import os
from math import floor, ceil, sqrt, log  

from ase.visualize import view



# Check if vasp path set correctly, if not, exit early
vasp_path = os.environ.get('EBROOTVASP') #EBROOTVASP, Niagara = SCINET_VASP_ROOT
if not isinstance(vasp_path, str):  
    exit()

num_cores = int(os.environ.get('SLURM_NPROCS'))
npar_setting = 2**ceil(log(floor(sqrt(num_cores)), 2)) # NPAR must be a divisor of the number of cores. Assumed to be a power of 2.

calc = vasp_calculator.Vasp(encut=500,
                        xc='PBE',
                        gga='PE',
                        ncore=4,
                        ivdw=12,
                        kpts  = (4,4,1),
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        sigma = 0.05,
                        nelm=400,
                        algo = 'fast',
                        ibrion=-1,    # -1 for no relaxation with vasp, 1 otherwise
                        ediffg=-0.01,  # forces
                        ediff=1e-5,  #energy conv.
                        prec='Accurate',
                        nsw=1, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=1 # 1 non-spin-polarized, #2 spin polarized
                        )


calc2 = vasp_calculator.Vasp(encut=400,
                        xc='PBE',
                        gga='PE',
                        ncore=4,
                        ivdw=12,
                        kpts  = (4,4,1),
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        sigma = 0.05,
                        nelm=400,
                        algo = 'fast',
                        ibrion=-1,    # -1 for no relaxation with vasp, 1 otherwise
                        ediffg=-0.01,  # forces
                        ediff=1e-5,  #energy conv.
                        prec='Accurate',
                        nsw=1, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=2 # 1 non-spin-polarized, #2 spin polarized
                        )

def print_atoms(atoms):

    print(atoms)
    print(atoms.get_chemical_formula())
    print(atoms.get_positions())


def relax(atoms, name):
    logfile = "{}_log.txt".format(name)
    traj_file = "{}_traj.traj".format(name)
    final_file = "{}.traj".format(name)

    atoms.set_calculator(calc)
    dyn = BFGS(atoms, logfile=logfile, trajectory=traj_file)
    dyn.run(fmax=0.01)
    write(final_file,atoms)


# Read the bulk structure
atoms = read("TaAgO3.poscar")
file_header = "bulk"

relax(atoms, "bulk")

# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 0, 0)  # Change this to (1,1,1) or any other plane

# Number of layers and vacuum thickness
layers = 1  # Adjust as needed
vacuum = 10.0  # Thickness of vacuum in Å

# Create the surface
slab = surface(atoms, miller_indices, layers, vacuum)


# Save
relax(atoms, "slab")

write('slab.traj', slab)

Atoms.calc = None
calc.set(ispin = 2)


# define adsorbates (OH and CO2 right now)
length = 0.96
OH = Atoms(['O', 'H'], charges=[-2, 1], positions=[[0, 0, 0], [0, 0, length]])

length = 1.162
CO2 = Atoms(['C', 'O', 'O'], positions=[[0, 0, 0], [0, 0, length / 3], [0, 0, 2*length/3]])
# Add OH adsorbate.

add_adsorbate(atoms, OH, 1.5)

relax(atoms, "OH")

# Save


# Close the trajectory file
# Relax


# Add CO2 adsorbate.


add_adsorbate(atoms, CO2, 1.5)

print_atoms(atoms)

relax(atoms, "CO2")

# relax