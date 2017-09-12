from ase import Atoms
from ase.calculators.nwchem import NWChem
from ase.db import connect
from ase.optimize import BFGS
from ase.io import read,write

import os

in_jsonfile  = "sode_series001_in.json"
out_jsonfile = "sode_series001_out.json"

# smiles = con1.get(abb_name="EC").smiles
num = 1
mol = read("%s@id=%d" % (in_jsonfile,num) )[0]

# neutral
calc_neu = NWChem(xc="B3LYP", charge=0, mult=1, iterations=100)
mol.set_calculator(calc_neu)

opt = BFGS(mol)
opt.run(fmax=0.05)

e_neu = mol.get_potential_energy()

"""
# cation
calc_cat = NWChem(xc="B3LYP", charge=+1, mult=2, iterations=100)
mol.set_calculator(calc_cat)

opt = BFGS(mol)
opt.run(fmax=0.05)

e_cat = mol.get_potential_energy()

# anion
calc_ani = NWChem(xc="B3LYP", charge=-1, mult=2, iterations=100)
mol.set_calculator(calc_ani)

opt = BFGS(mol)
opt.run(fmax=0.05)

e_ani = mol.get_potential_energy()

write("ec_opt.xyz",mol)

##########3

"""
con2 = connect(out_jsonfile)
con2.write(mol)

