from ase import Atoms
from ase.calculators.nwchem import NWChem

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write

import os, sys

ions = []
argvs = sys.argv

nions = len(argvs) - 2 # number of ion species

# get ion list
for i in range(0,nions):
	ions.append(str(argvs[i+1]))

print "ions", ions

num   = int(argvs[nions + 1])

orig_jsonfile     = "sode_orig.json"

# ---
xc    = "B3LYP"
basis = "3-21G"
fmax  =  0.10
# ---

db_ori   = connect(orig_jsonfile)

solv     = db_ori.get_atoms(num=num)
smiles   = db_ori.get(num=num).smiles
abb_name = db_ori.get(num=num).abb_name
wgt      = db_ori.get(num=num).molecular_weight
dens     = db_ori.get(num=num).density
bp       = db_ori.get(num=num).boiling_point
mp       = db_ori.get(num=num).melting_point
fp       = db_ori.get(num=num).flushing_point

#
# now loop over ions
#

for ion in ions:
	solv_ion_jsonfile = "sode_" + ion + ".json"

	label = "calc" + ion + str(num).zfill(4) + "/nwchem_low"

	# determine charge

	if ion == "Li":
		ion_charge = 1
	elif ion == "Na":
		ion_charge = 1
	elif ion == "Mg":
		ion_charge = 2
	else:
		print "only Li, Na, Mg is allowed" ; quit()

	db_ion = connect(solv_ion_jsonfile)

	print ion, "coordinated"

	replace_str1 = "([" + ion + "][O]="
	replace_str2 = "=[O][" + ion + "])"

	ion_smiles = "(" + smiles + ")"
	ion_smiles = ion_smiles.replace("(O=", replace_str1)
	ion_smiles = ion_smiles.replace("=O)", replace_str2)

	print "ion_smiles", ion_smiles

	babel_str = 'obabel ' + '-:"' + ion_smiles + '" -oxyz -O tmp.xyz -h --gen3D'
	os.system(babel_str)

	ion_solv  = read("tmp.xyz")

	ion_solv.calc = NWChem(label=label, xc=xc, basis=basis, charge=ion_charge, mult=1, 
						   iterations=200, mulliken=True) # cation
                       
	traj = ion + "_low_" + str(num).zfill(4) + ".traj"
	BFGS(ion_solv, trajectory=traj).run(fmax=fmax)

	db_ion.write(ion_solv, smiles=ion_smiles, abb_name=abb_name, num=num, 
					molecular_weight=wgt, density=dens,
					boiling_point=bp, melting_point=mp, flushing_point=fp,
				  	level="low"
			)

