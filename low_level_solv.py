from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.calculators.gaussian import Gaussian

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write

import os, sys

argvs = sys.argv
num = int(argvs[1])

calculator = "gaussian"
calculator = calculator.lower()

orig_jsonfile = "sode_orig.json"
solv_jsonfile = "tmp.json"

# ---
method = "B3LYP"
basis  = "6-31G"
fmax   = 0.10
# ---
label = "calc" + str(num).zfill(4) + "/low_solv"

db_ori = connect(orig_jsonfile)
db_sol = connect(solv_jsonfile)

solv     = db_ori.get_atoms(num=num)
smiles   = db_ori.get(num=num).smiles
abb_name = db_ori.get(num=num).abb_name
wgt      = db_ori.get(num=num).molecular_weight
dens     = db_ori.get(num=num).density
bp       = db_ori.get(num=num).boiling_point
mp       = db_ori.get(num=num).melting_point
fp       = db_ori.get(num=num).flushing_point

if "nw" in calculator:
	solv.calc = NWChem(label=label, xc=method, basis=basis, charge=0, mult=1, iterations=200,
			   mulliken=True)
elif "gau" in calculator:
	solv.calc = Gaussian(method=method, basis=basis, label=label, nprocshared=12)

traj = "solv_low_" + str(num).zfill(4) + ".traj"
BFGS(solv, trajectory=traj).run(fmax=fmax)

db_sol.write(solv, smiles=smiles, abb_name=abb_name, num=num, 
		   molecular_weight=wgt, density=dens,
		   boiling_point=bp, melting_point=mp, flushing_point=fp,
		   level="low"
	    )

