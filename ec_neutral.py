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

con1 = connect("in.json")

smiles = con1.get(name="ethylene carbonate").smiles

babel_str = 'obabel ' + '-:"' + smiles + '" -oxyz -O tmp.xyz -h --gen3D'

os.system(babel_str)

mol  = read("tmp.xyz")

# neutral
calc = NWChem(xc="B3LYP", basis=basis, charge=0, mult=1, iterations=100, \
              memory="2 gb", \
              mulliken=True)
mol.set_calculator(calc)

opt = BFGS(mol)
opt.run(fmax=fmax)

E_solv = mol.get_potential_energy()

mul_charge = calc.results["mul_charge"]

Ehomo = calc.get_homo_energy()
Elumo = calc.get_lumo_energy()

print "mulliken charge: ",mul_charge
print "homo energy: ", Ehomo
print "lumo energy: ", Elumo

### ----------- Li coordinated

print "Li-coordinated"
smiles  = con1.get(name="ethylene carbonate").smiles
smiles2 = "(" + smiles + ")"
li_smiles = smiles2.replace("(O=", "[Li](O=")

babel_str = 'obabel ' + '-:"' + li_smiles + '" -oxyz -O tmp2.xyz -h --gen3D'
os.system(babel_str)

mol  = read("tmp2.xyz")

# should be cation
calc = NWChem(xc="B3LYP", basis=basis, charge=1, mult=1, iterations=100, \
              mulliken=True)
mol.set_calculator(calc)

opt = BFGS(mol)
opt.run(fmax=fmax)

E_Li_solv = mol.get_potential_energy()

mul_charge = calc.results["mul_charge"]

e_homo = calc.get_homo_energy()
e_lumo = calc.get_lumo_energy()

print "mulliken charge: ",mul_charge
print "homo energy: ", e_homo
print "lumo energy: ", e_lumo

symbols = np.array( mol.get_chemical_symbols() )
Li_idx = np.where(symbols=="Li")[0]
O_list = np.where(symbols=="O")[0]
Li_O_dist = mol.get_distances(Li_idx, O_list)

R_Li_O = np.min(Li_O_dist)
coordinatingO = O_list[np.argmin(Li_O_dist)]

print "Li-O distance: ", R_Li_O

li = Atoms("Li", positions=[(0,0,0)])
calc = NWChem(xc="B3LYP", basis=basis, charge=1, mulliken=True,memory="total 2 gb")
li.set_calculator(calc)

E_Li = li.get_potential_energy()

Ecoord = E_Li_solv - ( E_Li + E_solv )

print "Li coordination energy (in eV)", Ecoord

