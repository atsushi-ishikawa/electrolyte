from ase import Atoms
# from ase.calculators.nwchem import NWChem
import mynwchem
from mynwchem import NWChem
from ase.db import connect
from ase.optimize import FIRE
from ase.io import read,write
import os, sys, shutil
import numpy as np

fmax = 0.20
xc = "B3LYP"
basis  = "3-21G"

db_read  = connect("ishi3.json")
db_write = connect("ishi3_new.json")

argvs = sys.argv
num = int(argvs[1])

# workdir
work = "/work/a_ishi/"
if not os.path.isdir(work):
	os.makedirs(work)

workdir = work + "calc" + str(num).zfill(4)
label   = workdir + "/low_solv"

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

charge = 0
if num == 106: # benzoate
	charge = -1
	
calc = NWChem(label=label, xc=xc, basis=basis, charge=charge, mult=1, iterations=100, \
              memory="8 gb", mulliken=True)
mol.set_calculator(calc)

opt = FIRE(mol)
opt.run(fmax=fmax)

db_write.write(mol, smiles=smiles, name=name, num=num,
               molecular_weight=wgt, density=dens,
               boiling_point=bp, melting_point=mp, flushing_point=fp, pubchemCID=pubchem
            )

shutil.rmtree(workdir)

