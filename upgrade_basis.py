from ase import Atoms
# from ase.calculators.nwchem import NWChem
from mynwchem import NWChem

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write

from tools import ABcoord

import os
import sys

#
# upgrade calculation of with specifiled jsonfile and num
#

argvs    = sys.argv
jsonfile = argvs[1]
num      = int(argvs[2])

# ---
xc     = "B3LYP"
basis  = "6-31G*"
fmax   =  0.10
memory = "total 2 gb"
# ---

label="calc" + str(num).zfill(4) + "/nwchem_upgrade"

db  = connect(jsonfile)

mol      = db.get_atoms(num=num)
smiles   = db.get(num=num).smiles
abb_name = db.get(num=num).abb_name
wgt      = db.get(num=num).molecular_weight
dens     = db.get(num=num).density
bp       = db.get(num=num).boiling_point
mp       = db.get(num=num).melting_point
fp       = db.get(num=num).flushing_point

# get charge and multiplicity from old file
charge = db.get(num=num).calculator_parameters["charge"]
mult   = db.get(num=num).calculator_parameters["mult"]

mol.calc = NWChem(label=label, xc=xc, basis=basis, charge=charge, mult=mult, iterations=200,
                  mulliken=True, memory=memory)

BFGS(mol).run(fmax=fmax)

db.write(mol, smiles=smiles, abb_name=abb_name, num=num, 
			  molecular_weight=wgt, density=dens,
			  boiling_point=bp, melting_point=mp, flushing_point=fp,
			  level="high"
			)

# delete old record
id = db.get(num=num, level="low").id
db.delete([id])

