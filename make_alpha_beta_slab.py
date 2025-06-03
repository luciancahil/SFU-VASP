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
from ase.build.tools import sort

import traceback






def fix_slab(old_slab, miller_index, remove_from_top, primitive_num_atoms, vacuum):
    slab = old_slab.repeat([1, 1, 1])
    epsilon = 4e-10

    if(miller_index == (1, 0, 0)):
        max_x = slab.cell[0][0]
        max_y = slab.cell[1][1]

        for atom in slab:
            if(atom.position[0] < epsilon):
                atom.position[0] = max_x
            
            if(atom.position[1] < epsilon):
                atom.position[1] = max_y
    

    slab = sort(slab, tags=slab.positions[:, 2])

    from_top_indicies = range(len(slab) - remove_from_top, len(slab))

    from_bot_indicies = range(primitive_num_atoms - remove_from_top)

    print("before: {}".format(len(slab)))
    if(remove_from_top != -1):
        if(remove_from_top != 0):
            del slab[[from_top_indicies] ]

        del slab[[from_bot_indicies] ]

    print("after: {}".format(len(slab)))

    slab.set_constraint(FixAtoms(indices=range(int(len(slab) / 2))))

    slab.center(vacuum, axis=2)

    return slab




# MAIN
alpha_remove = 1
beta_remove = 12

# Define the Miller indices (h, k, l) for the surface
miller_indices = (1, 0, 0)  # Change this to (1,1,1) or any other plane

layers = 3

name = "Bi2Ru2O7"

format = "poscar"

atoms = read("./POSCARS/{}.{}".format(name, format))



# Create the surface
dummy_slab = surface(atoms, miller_indices, 1, 0)
height = dummy_slab.cell[2][2]

primitive_num_atoms = len(atoms)

vacuum = 10

slab = surface(atoms, miller_indices, layers, vacuum=0)
slab_beta = surface(atoms, miller_indices, layers + 1, vacuum = 0)
slab_alpha = surface(atoms, miller_indices, layers + 1, vacuum = 0)

slab = fix_slab(slab, miller_indices, -1, primitive_num_atoms, vacuum)
slab_alpha = fix_slab(slab_alpha, miller_indices, alpha_remove, primitive_num_atoms, vacuum)
slab_beta = fix_slab(slab_beta, miller_indices, beta_remove, primitive_num_atoms, vacuum)

core_name = './slabs/{}_{}{}{}_slab'.format(name, miller_indices[0], miller_indices[1], miller_indices[2])
small_name = core_name + ".traj"
alpha_name = core_name + "_alpha.traj"
beta_name = core_name + "_beta.traj"


write(small_name, slab)
write(alpha_name, slab_alpha)
write(beta_name, slab_beta)

print("ase gui {}".format(small_name))
print("ase gui {}".format(alpha_name))
print("ase gui {}".format(beta_name))



"""
How to use.

1. Place the poscar of the primitive cell into the POSCARS folder.
2. Set "name" to the name of the file, minus the ".poscar" suffix. If you are using a different format, change the "format" variable
3. Set the miller indicies, amount of vacuum in Angstroms you want, and the number of layers you want.
4. Run once, and type the first "ase gui" command to see the slab generated.
5. Count how many atoms you want remvoed from the top for both the alpha and beta time. Set "alpha_remove" and "beta_remove" to those values (may often remove 0 from alpha).
6. Run again, and check your final structures to see if you're happy with them


"""