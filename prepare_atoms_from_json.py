from ase import Atoms
# from ase.calculators.nwchem import NWChem
from ase.calculators.gaussian import Gaussian
import mynwchem
from mynwchem import NWChem
from ase.db import connect
from ase.optimize import FIRE
from ase.io import read,write
from tools import delete_num_from_json
import os, sys, shutil
import numpy as np

fmax  = 0.20
xc    = "B3LYP"
basis = "3-21G"
calculator = "gaussian"
opt   = "gediis"

orig_jsonfile = "ishi3.json"
solv_jsonfile = "ishi3_new.json"

db_read  = connect(orig_jsonfile)
db_write = connect(solv_jsonfile)

argvs = sys.argv
num = int(argvs[1])

# workdir
work    = "/work/a_ishi/"
workdir = work + "calc" + str(num).zfill(4)
#label   = workdir + "/low_solv"
label = workdir + "/prep"

if not os.path.isdir(label):
	os.makedirs(label)

smiles  = db_read.get(num=num).smiles
name    = db_read.get(num=num).name
wgt     = db_read.get(num=num).molecular_weight
dens    = db_read.get(num=num).density
bp      = db_read.get(num=num).boiling_point
mp      = db_read.get(num=num).melting_point
fp      = db_read.get(num=num).flushing_point
pubchem = db_read.get(num=num).pubchemCID

tmpxyz = label + "/tmp.xyz"
babel_str = 'obabel ' + '-:"' + smiles + '" -oxyz -O ' + tmpxyz + ' -h --gen3D'
os.system(babel_str)

mol = read(tmpxyz)

charge = 0
if num == 106: # benzoate
	charge = -1
	
if "gau" in calculator:
	mol.calc = Gaussian(label=label, method=xc, basis=basis, charge=charge, mult=1, force=None, opt=opt)
	mol.get_potential_energy()
else:
	mol.calc = NWChem(label=label, xc=xc, basis=basis, charge=charge, mult=1, iterations=100, memory="8 gb", mulliken=True)
	FIRE(mol).run(fmax=fmax)

delete_num_from_json(num, solv_jsonfile)
db_write.write(mol, smiles=smiles, name=name, num=num,
               molecular_weight=wgt, density=dens,
               boiling_point=bp, melting_point=mp, flushing_point=fp, pubchemCID=pubchem )

shutil.rmtree(workdir)

