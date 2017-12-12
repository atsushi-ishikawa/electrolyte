from ase import Atoms
# from ase.calculators.nwchem import NWChem
import mynwchem
from mynwchem import NWChem

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read,write

import os
import numpy as np

fmax = 0.05
basis  = "sto-3g"

db_read  = connect("ishi3.json")
db_write = connect("ishi3_new.json")

num = 1

smiles  = db_read.get(num=num).smiles
name    = db_read.get(num=num).name
wgt     = db_read.get(num=num).molecular_weight
dens    = db_read.get(num=num).density
bp      = db_read.get(num=num).boiling_point
mp      = db_read.get(num=num).melting_point
fp      = db_read.get(num=num).flushing_point
pubchem = db_read.get(num=num).pubchemCID

babel_str = 'obabel ' + '-:"' + smiles + '" -oxyz -O tmp.xyz -h --gen3D'

os.system(babel_str)

mol  = read("tmp.xyz")

calc = NWChem(xc="B3LYP", basis=basis, charge=0, mult=1, iterations=100, \
              memory="2 gb", \
              mulliken=True)
mol.set_calculator(calc)

opt = BFGS(mol)
opt.run(fmax=fmax)

db_write.write(mol, smiles=smiles, name=name, num=num,
               molecular_weight=wgt, density=dens,
               boiling_point=bp, melting_point=mp, flushing_point=fp, pubchemCID=pubchem
            )

