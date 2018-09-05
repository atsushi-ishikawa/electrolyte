from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.calculators.gaussian import Gaussian
from tools import delete_num_from_json
from ase.db import connect
from ase.optimize import FIRE
from ase.io import read, write
import os, sys, shutil

argvs = sys.argv
num = int(argvs[1])

calculator = "gaussian"
calculator = calculator.lower()

#solv_jsonfile = "electrolye_2017Aug.json"
solv_jsonfile = "ishi3_new.json" # read from & write to this file

method = "B3LYP"
basis  = "3-21G"
fmax   = 0.10

# workdir
work = "/work/a_ishi/"
if not os.path.isdir(work):
	os.makedirs(work)

workdir = work + "calc" + str(num).zfill(4)
label   = workdir + "/low_solv"

db_solv = connect(solv_jsonfile)

solv   = db_solv.get_atoms(num=num)
smiles = db_solv.get(num=num).smiles
name   = db_solv.get(num=num).name
wgt    = db_solv.get(num=num).molecular_weight
dens   = db_solv.get(num=num).density
bp     = db_solv.get(num=num).boiling_point
mp     = db_solv.get(num=num).melting_point
fp     = db_solv.get(num=num).flushing_point
try:
	pubchem = db_read.get(num=num).pubchemCID
except:
	pubchem = ""


charge = 0
if num == 106:
	charge = -1

if "gau" in calculator:
	solv.calc = Gaussian(method=xc, basis=basis, label=label, 
						 charge=charge, mult=1, nprocshared=12, population="full")

elif "nw" in calculator:
	solv.calc = NWChem(label=label, xc=method, basis=basis, charge=charge, mult=1, iterations=200,
			   mulliken=True)

traj = "solv_low_" + str(num).zfill(4) + ".traj"
FIRE(solv, trajectory=traj).run(fmax=fmax)

delete_num_from_json(num, solv_jsonfile)

db_solv.write(solv, smiles=smiles, name=name, num=num, 
			  molecular_weight=wgt, density=dens,
 			  boiling_point=bp, melting_point=mp, flushing_point=fp, pubchemCID=pubchem
 			 )

shutil.rmtree(workdir)

