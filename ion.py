from ase import Atoms
# from ase.calculators.nwchem import NWChem
import mynwchem
from mynwchem import NWChem
# from ase.calculators.gaussian import Gaussian

from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write

from ase.visualize import view

from tools import delete_num_from_json

import os, sys

ions = []
argvs = sys.argv

calculator = "nwchem"
calculator = calculator.lower()

nions = len(argvs) - 2 # number of ion species

# get ion list
for i in range(0, nions):
	ions.append(str(argvs[i+1]))

print "ions", ions

num = int(argvs[nions + 1])

# solv_jsonfile = "electrolye_2017Aug.json"
solv_jsonfile = "ishi3_new.json"

# ---
xc    = "B3LYP"
#basis = "6-31G*"
basis = "DZP"
fmax  =  0.10
memory = "total 8 gb"
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
try:
	pubchem = db_solv.get(num=num).pubchemCID
except:
	pubchem = ""

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

		if ion in ["Li", "Na", "K", "Rb", "Cs"]:
			ion_charge = 1
		elif ion == "Mg":
			ion_charge = 2
		else:
			print "only AlcaliMetal and Mg is allowed" ; quit()

		print "smiles", smiles
		if "(O=" in smiles or "=O)" in smiles:
			print "case 1"
			ion_smiles = "" + smiles + ""
			replace_str1 = "(([" + ion + "])O="
			replace_str2 = "=O([" + ion + "]))"
			ion_smiles = ion_smiles.replace("(O=", replace_str1,1)
			ion_smiles = ion_smiles.replace("=O)", replace_str2,1)
		elif "O=" in smiles or "=O" in smiles:
			print "case 2"
			ion_smiles = "" + smiles + ""
			replace_str1 = "([" + ion + "])O="
			replace_str2 = "=O([" + ion + "])"
			ion_smiles = ion_smiles.replace("O=", replace_str1,1)
			ion_smiles = ion_smiles.replace("=O", replace_str2,1)
		elif "[O]" in smiles:
 			print "case 3c"
 			ion_smiles = "" + smiles + ""
 			replace_str1 = "[O]([" + ion + "])"
 			ion_smiles = ion_smiles.replace("[O]", replace_str1,1)
		elif "O1" in smiles:
			if smiles.startswith("O1"):
				print "case 3-A"
				ion_smiles = smiles
				replace_str1 = "[" + ion + "]O1"
		 		#replace_str1 = "O1"
				ion_smiles = ion_smiles.replace("O1", replace_str1,1)
			else:
				print "case 3-a"
				ion_smiles = smiles
				replace_str1 = "O1([" + ion + "])"
		 		#replace_str1 = "O1"
				ion_smiles = ion_smiles.replace("O1", replace_str1,1)
		elif "o1" in smiles:
			print "case 3-b"
			ion_smiles = smiles
			#replace_str1 = "[" + ion + "](o1)"
			replace_str1 = "o1([" + ion + "])"
			ion_smiles = ion_smiles.replace("o1", replace_str1,1)
		elif "O" in smiles:
 			print "case 3d"
 			ion_smiles = "" + smiles + ""
 			#replace_str1 = "[" + ion + "]O"
 			replace_str1 = "O([" + ion + "])"
 			ion_smiles = ion_smiles.replace("O", replace_str1,1)
		elif "C#N" in smiles:
			print "case 4"
			ion_smiles = smiles
			replace_str1 = "C#N([" + ion + "])"
			ion_smiles = ion_smiles.replace("C#N", replace_str1,1)
		elif "n1" in smiles:
			print "case 5-a"
			ion_smiles = smiles
			replace_str1 = "n1([" + ion + "])"
			ion_smiles = ion_smiles.replace("n1", replace_str1,1)
		elif "n2" in smiles:
			print "case 5-b"
			ion_smiles = smiles
			replace_str1 = "n2([" + ion + "])"
			ion_smiles = ion_smiles.replace("n2", replace_str1,1)
		elif "S" in smiles:
			print "case 6-a"
			ion_smiles = smiles
			replace_str1 = "S([" + ion + "])"
			ion_smiles = ion_smiles.replace("S", replace_str1,1)
		elif "s1" in smiles:
			print "case 6-b"
			ion_smiles = smiles
			replace_str1 = "s1([" + ion + "])"
			ion_smiles = ion_smiles.replace("s1", replace_str1,1)
		else:
			print "case 10"
			ion_smiles = "(" + smiles + ")"
			ion_smiles = ion_smiles + "[" + ion + "]"

#		print ion, "coordinated"
#
#		replace_str1 = "([" + ion + "][O]="
#		replace_str2 = "=[O][" + ion + "])"
#
#		ion_smiles = "(" + smiles + ")"
#		ion_smiles = ion_smiles.replace("(O=", replace_str1, 1) # replace just one
#		ion_smiles = ion_smiles.replace("=O)", replace_str2, 1) # 

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

	if "nw" in calculator:
	 	ion_solv.calc = NWChem(label=label, xc=xc, basis=basis, charge=ion_charge, mult=1, 
 				       iterations=200, mulliken=True, memory=memory) # cation
	elif "gau" in calculator:
 		ion_solv.calc = Gaussian(method=xc, basis=basis, label=label, 
					 charge=ion_charge, mult=1, nprocshared=12)
                        
 	traj = ion + "_low_" + str(num).zfill(4) + ".traj"
 	BFGS(ion_solv, trajectory=traj).run(fmax=fmax)
 
 	db_ion.write(ion_solv, smiles=ion_smiles, name=name, num=num, 
				molecular_weight=wgt, density=dens,
 				boiling_point=bp, melting_point=mp, flushing_point=fp, pubchemCID=pubchem,
 			  	level="low"
 			)

