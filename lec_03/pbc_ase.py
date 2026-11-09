import numpy as np
from ase import Atoms,units
from ase.io import write
from ase.md.velocitydistribution import thermalize_momenta,Stationary
from ase.md.verlet import VelocityVerlet
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS
import matplotlib.pyplot as plt

box_l = 10
box = [box_l,box_l,box_l]
atoms = Atoms(symbols = ['Ar','Ar'], positions =[[0,0,0],[0,0,3.2]],cell = box ,pbc = True)
thermalize_momenta(atoms,temperature_K = 300)
atoms.calc = LennardJones(sigma = 3.4,epsilon =0.014,rc = 5.5)
opt = BFGS(atoms)
opt.run(fmax=0.02)
atoms.wrap()
md = VelocityVerlet(atoms,timestep = 1*units.fs,logfile ='-',trajectory = 'traj.md')
md.run(10000)