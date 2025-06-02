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

name = "TaAgO3"

type = "poscar"

atoms = read("./POSCARS/{}.{}".format(name, type))



# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 0, 0)  # Change this to (1,1,1) or any other plane

# TODO: Fix cell, duplicate in Y direction.
# Number of layers and vacuum thickness
layers = 6  # Adjust as needed
vacuum = 10.0  # Thickness of vacuum in Å
# Create the surface

slab = surface(atoms, miller_indices, layers, vacuum)




slab = slab.repeat([1, 2, 1])
epsilon = 4e-10

max_x = slab.cell[0][0]
max_y = slab.cell[1][1]

for atom in slab:
    if(atom.position[0] < epsilon):
        atom.position[0] = max_x
    
    if(atom.position[1] < epsilon):
        atom.position[1] = max_y

z_positions = slab.get_positions()[:, 2]
threshold = sorted(z_positions)[int(len(z_positions) * 1 / 2)]  # freeze bottom 2 of 4
frozen_indices = [i for i, z in enumerate(z_positions) if z < threshold]
slab.set_constraint(FixAtoms(indices=frozen_indices))

slab.pbc = [True, True, True]

write('{}_{}{}{}_slab.traj'.format(name, miller_indices[0], miller_indices[1], miller_indices[2]), slab)
write('{}_{}{}{}_slab_beta.traj'.format(name, miller_indices[0], miller_indices[1], miller_indices[2]), slab)

