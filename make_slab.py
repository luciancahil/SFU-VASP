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


def fix_slab(old_slab, miller_index):
    slab = old_slab.repeat([1, 2, 1])
    epsilon = 4e-10

    if(miller_index == (1, 0, 0)):
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

name = "Bi2Ru2O7"

type = "poscar"

atoms = read("./POSCARS/{}.{}".format(name, type))


# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 1, 1)  # Change this to (1,1,1) or any other plane

# TODO: Fix cell, duplicate in Y direction.
# Number of layers and vacuum thickness
layers = 4  # Adjust as needed
vacuum = 10.0  # Thickness of vacuum in Å
# Create the surface
dummy_slab = surface(atoms, miller_indices, 1, 0)
height = dummy_slab.cell[2][2]

slab = surface(atoms, miller_indices, layers, vacuum=0)
slab_beta = surface(atoms, miller_indices, layers + 1, vacuum = 0)

# filter out the top half of the top cell and bottom half of the bottom cell

index_height_pair = [(a.index, a.position[2].item()) for a in slab_beta]

index_height_pair = sorted(index_height_pair, key= lambda x:x[1])




slab = fix_slab(slab, miller_indices)

slab_beta = fix_slab(slab_beta, miller_indices)







write('{}_{}{}{}_slab.traj'.format(name, miller_indices[0], miller_indices[1], miller_indices[2]), slab)
write('{}_{}{}{}_slab_extra.traj'.format(name, miller_indices[0], miller_indices[1], miller_indices[2]), slab_beta)
print(slab.cell)
print("ase gui {}_{}{}{}_slab.traj".format(name, miller_indices[0], miller_indices[1], miller_indices[2]))
print("ase gui {}_{}{}{}_slab_extra.traj".format(name, miller_indices[0], miller_indices[1], miller_indices[2]))
