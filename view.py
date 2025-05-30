#! /usr/bin/env python

from ase.io import read, Trajectory, write
from ase.build import fcc111, root_surface,add_adsorbate, add_vacuum
from ase import Atoms
import ase.calculators.vasp as vasp_calculator
from ase.optimize import BFGS
from ase.build import surface
import os
from math import floor, ceil, sqrt, log  
from ase.constraints import FixAtoms

from ase.visualize import view



# Read the bulk structure
atoms = read("TaAgO3.poscar")
file_header = "bulk"


# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 0, 0)  # Change this to (1,1,1) or any other plane

# Number of layers and vacuum thickness
layers = 4  # Adjust as needed
vacuum = 15.0  # Thickness of vacuum in Å

# Create the surface
slab = surface(atoms, miller_indices, layers, vacuum)



# Freeze the bottom two layers
z_positions = slab.get_positions()[:, 2]
threshold = sorted(z_positions)[int(len(z_positions) * 1 / 2)]  # freeze bottom 2 of 6
frozen_indices = [i for i, z in enumerate(z_positions) if z < threshold]
slab.set_constraint(FixAtoms(indices=frozen_indices))

slab.pbc = [True, True, True]

breakpoint()




# Save

write('slab.traj', slab)

Atoms.calc = None


# define adsorbates (OH and CO2 right now)
length = 0.96
OH = Atoms(['O', 'H'], charges=[-2, 1], positions=[[0, 0, 0], [0, 0, length]])

length = 1.162
CO2 = Atoms(['C', 'O', 'O'], positions=[[0, 0, 0], [0, 0, length / 3], [0, 0, 2*length/3]])
# Add OH adsorbate.

add_adsorbate(slab, OH, 1.5)


write('slabOH.traj', slab)
breakpoint()

# Save


# Close the trajectory file
# Relax


# Add CO2 adsorbate.


add_adsorbate(atoms, CO2, 1.5)
write('slabCO2.traj', slab)
view(atoms)



# relax


