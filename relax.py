#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import ase
#from ase.constraints import FixAtoms
#from ase.structure import bulk
#from ase.lattice.spacegroup import crystal
from ase.io import *
import sys 
import os
#if not (os.path.exists('vdw_kernel.bindat')):
#                os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
import numpy as np

#atoms=read('in.cif')
#atoms=read('CONTCAR')

try:
  atoms=read('slab_B.traj')
except:
  atoms=read('in.xyz')
  atoms.center(vacuum=9)
  atoms.pbc=True
try:
  submitdir=sys.argv[1]
except:
  submitdir=''
if submitdir != '':
  submitdir += '/'

calc = vasp_calculator.Vasp(encut=400,
                        xc='PBE',
                        gga='PE',
                        ncore=4,
                        ivdw=11,
                        kpts  = (4,4,1),
                        gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        ismear=0,
                        sigma = 0.1,
                        nelm=400,
                        algo = 'fast',
                        ibrion=-1,    # -1 for no relaxation with vasp, 1 otherwise
                        ediffg=-0.01,  # forces
                        ediff=1e-8,  #energy conv.
                        prec='Accurate',
                        nsw=1, # don't use the VASP internal relaxation, only use ASE
                        lreal='Auto',
                        ispin=1)
atoms.set_calculator(calc)
trajfile='final.traj'
dyn = BFGS(atoms, logfile=submitdir+'relax.log', trajectory=submitdir+'relax.traj')
dyn.run(fmax=0.01)
write(trajfile,atoms)
relaxed_atoms = read(trajfile)
e=relaxed_atoms.get_potential_energy()
print('final energy',e)
f=open('final.e','w')
f.write(str(e)+'\n')
f.close()
