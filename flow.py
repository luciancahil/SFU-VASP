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
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--crystal', type=str, default=None, help='The path to crystal we want to relax.', required=True)
parser.add_argument('--settings', type=str, default=None, help='The path to the settings folder.')

args = parser.parse_args()


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


# TODO: Process. Get files.

print("Warning: This file should never be run from the main folder. It should always be run from a child folder of the main folder.")
settings = open("../../processing/{}".format(args.settings), mode='r')

magmoms = None

for line in settings:
    parts = line.strip().split(",")
    key = parts[0]
    val = parts[1]


    if(key == "kpts"):
        val = int(val)
        calc_settings[key] = (val, val, 1)
    elif(key == "magmom"):
        if(val == -1):
            raise("The magmom is not configure properly. Please see the README for how to configure")

        magmoms = dict()

        pairs = val.split(" ")
        for pair in pairs:

            magmoms[pair.split(":")[0]] = float(pair.split(":")[1])
    else:
        if(val == "True" or val == "False"):
            val = bool(val)
        else:
            try:
                val = float(val)
            except:
                pass
            
        print(val)
        calc_settings[key] = val

calc = vasp_calculator.Vasp(**calc_settings)

atoms = read("../../POSCARS/{}".format(args.crystal))
file_header = "bulk"

if(magmoms != None):
    atoms.set_initial_magnetic_moments([magmoms[atom.symbol] for atom in atoms])

atoms.pbc = True

# Check if vasp path set correctly, if not, exit early
vasp_path = os.environ.get('EBROOTVASP') #EBROOTVASP, Niagara = SCINET_VASP_ROOT






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



#relax(atoms, "bulk")

try:
    relax(atoms, "slab")
except(ValueError):
    print(traceback.format_exc())
    calc.set(sigma = 0.2)
    relax(atoms, "slab")

write('slab.traj', atoms)




"""
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
"""