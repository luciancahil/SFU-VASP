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
from ase.build import cut
import traceback


name = "TaAgO3"

type = "poscar"

atoms = read("./POSCARS/{}.{}".format(name, type))


slab = cut(atoms,
           a=(1, -1, 0),
           b=(0, 1, -1),
           c=(1, 1, 1),
           nlayers=8)

#slab = slab.repeat([1, 2, 2])

write("slab_cut.traj", slab)