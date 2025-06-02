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


def fix_slab(old_slab):
    slab = old_slab.repeat([1, 2, 1])
    epsilon = 4e-10

    max_x = slab.cell[0][0]
    max_y = slab.cell[1][1]

    for atom in slab:
        if(atom.position[0] < epsilon):
            atom.position[0] = max_x
        
        if(atom.position[1] < epsilon):
            atom.position[1] = max_y

    z_positions = slab.get_positions()[:, 2]
    sorted_z = sorted(z_positions)
    halfway_atom = int(len(slab) / 2)
    
    # take the average of the 2 middle ones as a cutoff
    threshold = (sorted_z[halfway_atom] + sorted_z[halfway_atom - 1]) / 2
    
    frozen_indices = [i for i, z in enumerate(z_positions) if z <= threshold]

    slab.set_constraint(FixAtoms(indices=frozen_indices))

    slab.pbc = [True, True, True]

    return slab

name = "TaAgO3"

type = "poscar"

atoms = read("./POSCARS/{}.{}".format(name, type))

height = atoms.cell[2][2]

print(atoms.cell)
# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 0, 0)  # Change this to (1,1,1) or any other plane

# TODO: Fix cell, duplicate in Y direction.
# Number of layers and vacuum thickness
layers = 4  # Adjust as needed
vacuum = 10.0  # Thickness of vacuum in Å
# Create the surface

slab = surface(atoms, miller_indices, layers, vacuum)
slab_beta = surface(atoms, miller_indices, layers + 1, vacuum - height / 2)


slab = fix_slab(slab)

slab_beta = fix_slab(slab_beta)


# filter out the top half of the top cell and bottom half of the bottom cell
z_pos =  (slab_beta.get_positions()[:, 2])
lowest_pos = min(z_pos) + height / 4
highest_pos = max(z_pos) - height / 4

print(lowest_pos)
print(highest_pos)

([(atom.index, atom.position[2]) for atom in slab_beta if (atom.position[2] > highest_pos or atom.position[2] < lowest_pos)])

del slab_beta[[atom.index for atom in slab_beta if (atom.position[2] > highest_pos or atom.position[2] < lowest_pos)]]




write('{}_{}{}{}_slab.traj'.format(name, miller_indices[0], miller_indices[1], miller_indices[2]), slab)
write('{}_{}{}{}_slab_beta.traj'.format(name, miller_indices[0], miller_indices[1], miller_indices[2]), slab_beta)

print("ase gui {}_{}{}{}_slab.traj".format(name, miller_indices[0], miller_indices[1], miller_indices[2]))
print("ase gui {}_{}{}{}_slab_beta.traj".format(name, miller_indices[0], miller_indices[1], miller_indices[2]))
