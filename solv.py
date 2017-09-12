from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.calculators.gaussian import Gaussian
from tools import delete_num_from_json

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write

import os, sys

argvs = sys.argv
num = int(argvs[1])

calculator = "nwchem"
calculator = calculator.lower()

solv_jsonfile = "electrolye_2017Aug.json"

# ---
method = "B3LYP"
basis  = "6-31G"
fmax   = 0.10
# ---
label = "calc" + str(num).zfill(4) + "/low_solv"

db_solv = connect(solv_jsonfile)

solv   = db_solv.get_atoms(num=num)
smiles = db_solv.get(num=num).smiles
name   = db_solv.get(num=num).name
wgt    = db_solv.get(num=num).molecular_weight
dens   = db_solv.get(num=num).density
bp     = db_solv.get(num=num).boiling_point
mp     = db_solv.get(num=num).melting_point
fp     = db_solv.get(num=num).flushing_point

if "nw" in calculator:
	solv.calc = NWChem(label=label, xc=method, basis=basis, charge=0, mult=1, iterations=200,
			   mulliken=True)
elif "gau" in calculator:
	solv.calc = Gaussian(method=method, basis=basis, label=label, nprocshared=12)

traj = "solv_low_" + str(num).zfill(4) + ".traj"
BFGS(solv, trajectory=traj).run(fmax=fmax)

###
delete_num_from_json(num, solv_jsonfile)

db_solv.write(solv, smiles=smiles, name=name, num=num, 
		molecular_weight=wgt, density=dens,
		boiling_point=bp, melting_point=mp, flushing_point=fp
	    )

