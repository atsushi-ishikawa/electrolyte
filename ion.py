from ase import Atoms
from ase.calculators.nwchem import NWChem

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write

from tools import delete_num_from_json

import os, sys

ions = []
argvs = sys.argv

nions = len(argvs) - 2 # number of ion species

# get ion list
for i in range(0,nions):
	ions.append(str(argvs[i+1]))

print "ions", ions

num = int(argvs[nions + 1])

solv_jsonfile     = "electrolye_2017Aug.json"

# ---
xc    = "B3LYP"
basis = "6-31G*"
fmax  =  0.10
# ---

db_solv = connect(solv_jsonfile)

solv   = db_solv.get_atoms(num=num)
smiles = db_solv.get(num=num).smiles
name   = db_solv.get(num=num).name
wgt    = db_solv.get(num=num).molecular_weight
dens   = db_solv.get(num=num).density
bp     = db_solv.get(num=num).boiling_point
mp     = db_solv.get(num=num).melting_point
fp     = db_solv.get(num=num).flushing_point

#
# now loop over ions
#

for ion in ions:
	solv_ion_jsonfile = solv_jsonfile.split(".")[0] + "_" + ion + ".json"
	db_ion = connect(solv_ion_jsonfile)

	try:
		ion_solv = db_ion.get_atoms(num=num)
	except:
		print "newly formed here"
		# determine charge

		if ion == "Li":
			ion_charge = 1
		elif ion == "Na":
			ion_charge = 1
		elif ion == "Mg":
			ion_charge = 2
		else:
			print "only Li, Na, Mg is allowed" ; quit()


		print ion, "coordinated"

		replace_str1 = "([" + ion + "][O]="
		replace_str2 = "=[O][" + ion + "])"

		ion_smiles = "(" + smiles + ")"
		ion_smiles = ion_smiles.replace("(O=", replace_str1, 1) # replace just one
		ion_smiles = ion_smiles.replace("=O)", replace_str2, 1) # 

		print "ion_smiles", ion_smiles

		babel_str = 'obabel ' + '-:"' + ion_smiles + '" -oxyz -O tmp.xyz -h --gen3D'
		os.system(babel_str)

		ion_solv  = read("tmp.xyz")
	else:
		ion_charge = db_ion.get(num=num).calculator_parameters["charge"]
		ion_smiles = db_ion.get(num=num).smiles
		print "found in database"
		delete_num_from_json(num, solv_ion_jsonfile)

	label = "calc" + ion + str(num).zfill(4) + "/nwchem_low"

	ion_solv.calc = NWChem(label=label, xc=xc, basis=basis, charge=ion_charge, mult=1, 
						   iterations=200, mulliken=True) # cation
                       
	traj = ion + "_low_" + str(num).zfill(4) + ".traj"
	BFGS(ion_solv, trajectory=traj).run(fmax=fmax)

	db_ion.write(ion_solv, smiles=ion_smiles, name=name, num=num, 
					molecular_weight=wgt, density=dens,
					boiling_point=bp, melting_point=mp, flushing_point=fp,
				  	level="low"
			)

