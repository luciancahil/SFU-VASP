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
from ase.constraints import FixAtoms

import traceback
import numpy as np


num_cores = int(os.environ.get('SLURM_NPROCS', 1))
npar_setting = 2**ceil(log(floor(sqrt(num_cores)), 2)) # NPAR must be a divisor of the number of cores. Assumed to be a power of 2.


# Default settings
calc_settings = {
    'encut':500,
    'xc':'PBE',
    'gga':'PE',
    'ncore':4,
    'ivdw':12,
    'kpts'  : (4,4,1),
    'gamma' : True, # Gamma-centered (defaults to Monkhorst-Pack)
    'ismear':0,
    'sigma' : 0.05,
    'nelm':400,
    'algo' : 'fast',
    'ibrion':-1,    # -1 for no relaxation with vasp, 1 otherwise
    'ediffg':-0.01,  # forces
    'ediff':1e-5,  #energy conv.
    'prec':'Accurate',
    'nsw':1, # don't use the VASP internal relaxation, only use ASE
    'lreal':'Auto',
    'ispin':-1 # 1 non-spin-polarized, #2 spin polarized
}

print("Hi!")


# TODO: Process. Get files.

settings = open("processing/settings.csv", mode='r')


for line in settings:
    parts = line.strip().split(",")
    key = parts[0]
    try:
        val = float(parts[1])
    except:
        pass


    if(key == "kpoints"):
        val = int(val)
        calc_settings[key] = (val, val, 1)
    else:
        calc_settings[key] = val

calc = vasp_calculator.Vasp(**calc_settings)


# Check if vasp path set correctly, if not, exit early
vasp_path = os.environ.get('EBROOTVASP') #EBROOTVASP, Niagara = SCINET_VASP_ROOT
print("Vasp_path: {}".format(vasp_path))
if not isinstance(vasp_path, str):  
    exit()





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

#relax(atoms, "bulk")

# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 0, 0)  # Change this to (1,1,1) or any other plane

# TODO: Fix cell, duplicate in Y direction.
# Number of layers and vacuum thickness
layers = 4  # Adjust as needed
vacuum = 15.0  # Thickness of vacuum in Å
breakpoint()
# Create the surface
slab = surface(atoms, miller_indices, layers, vacuum)


write('slab.traj', slab)

1/0
slab.pbc = [True, True, True]


# Freeze the bottom two layers
z_positions = slab.get_positions()[:, 2]
threshold = sorted(z_positions)[int(len(z_positions) * 1 / 2)]  # freeze bottom 2 of 4
frozen_indices = [i for i, z in enumerate(z_positions) if z < threshold]
slab.set_constraint(FixAtoms(indices=frozen_indices))

try:
    relax(slab, "slab")
except(ValueError):
    print(traceback.format_exc())
    calc.set(sigma = 0.2)
    relax(slab, "slab")

write('slab.traj', slab)

Atoms.calc = None
calc.set(ispin = 2)


# define adsorbates (OH and CO2 right now)
length = 0.96
OH = Atoms(['O', 'H'], charges=[-2, 1], positions=[[0, 0, 0], [0, 0, length]])

length = 1.162
CO2 = Atoms(['O', 'C', 'O'], positions=[[0, 0, 0], [0, 0, length / 3], [0, 0, 2*length/3]])
# Add OH adsorbate.

add_adsorbate(slab, OH, 1.5)

relax(atoms, "OH")

# Save


# Close the trajectory file
# Relax


# Add CO2 adsorbate.


add_adsorbate(slab, CO2, 1.5)

print_atoms(atoms)

relax(atoms, "CO2")

# relax
