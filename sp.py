from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.db import connect
from ase.io import read,write
from ase.optimize import BFGS
import os, sys

argvs   = sys.argv
inxyz   = argvs[1]
outjson = "singlepoint.json"
db = connect(outjson)

# ---
xc     = "B3LYP"
basis  = "3-21G"
fmax   =  0.05
memory = "total 8 gb"
charge = 1
mult   = 1
# ---

mol = read(inxyz)
mol.calc = NWChem(xc=xc, basis=basis, charge=charge, mult=mult, iterations=200, memory=memory)
BFGS(mol, trajectory="tmp.traj").run(fmax=fmax)

db.write(mol)

