from ase import Atoms

from ase.db import connect
from ase.optimize import BFGS, FIRE
from ase.io import read, write
from ase.visualize import view

import os, sys

argvs  = sys.argv
dbfile = argvs[1]
num = int(argvs[2])

db  = connect(dbfile)
mol = db.get_atoms(num=num)

view(mol)

