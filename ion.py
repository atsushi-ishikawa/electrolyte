from ase import Atoms
# from ase.calculators.nwchem import NWChem
from mynwchem import NWChem
from mygaussian import *
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from ase.optimize import FIRE
from ase.io import read, write
from ase.visualize import view
from tools import delete_num_from_json
import os, sys, shutil

ions = []
argvs = sys.argv

calculator = "gaussian"
calculator = calculator.lower()

nions = len(argvs) - 2 # number of ion species

# get ion list
for i in range(0, nions):
	ions.append(str(argvs[i+1]))

num = int(argvs[nions + 1])

# solv_jsonfile = "electrolye_2017Aug.json"
solv_jsonfile = "ishi3_new.json"

# ---
xc     = "M06"
basis  = "lanl2dz"
pseudo = "lanl2"
#basis = "3-21G"
# ---
fmax    = 0.08
steps   = 200
#memory  = "total 8 gb"
opt     = "newton, maxcycles=200, loose"
scf     = "(xqc,maxconventional=200)"
grid    = "ultrafine"
ioplist = ["1/18=40"] # do optimization in Cartesian coordinate to avoid error
nprocs  = 12
mem     = "8GB"
# ---

# workdir
work    = "/work/a_ishi/"
workdir = work + "calc" + str(num).zfill(4)

if not os.path.isdir(workdir):
	os.makedirs(workdir)

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
		# look the database for the same solvent
		ion_solv = db_ion.get_atoms(num=num)
	except:
		print "newly formed here"
		# not found in the database -- newly form
		try:
			#
			# make coordination by replacing Li-complex
			#
			solv_li_jsonfile = solv_jsonfile.split(".")[0] + "_" + "Li" + ".json"
			db_li   = connect(solv_li_jsonfile)
			li_solv = db_li.get_atoms(num=num)
			symbols = li_solv.get_chemical_symbols()
			for i,j in enumerate(symbols):
				symbols[i] = j.replace("Li",ion)
			print "found Li then replace"
			li_solv.set_chemical_symbols(symbols)
			ion_solv = li_solv
		except:
			#
			# no previons information -- make from SMILES string
			#
			print "make from smiles = ", smiles
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
 				print "case 3"
 				ion_smiles = "" + smiles + ""
 				replace_str1 = "[O]([" + ion + "])"
 				ion_smiles = ion_smiles.replace("[O]", replace_str1,1)
			elif "O1" in smiles:
				if smiles.startswith("O1"):
					print "case 4A"
					ion_smiles = smiles
					replace_str1 = "[" + ion + "]O1"
			 		#replace_str1 = "O1"
					ion_smiles = ion_smiles.replace("O1", replace_str1,1)
				else:
					print "case 4B"
					ion_smiles = smiles
					replace_str1 = "O1([" + ion + "])"
			 		#replace_str1 = "O1"
					ion_smiles = ion_smiles.replace("O1", replace_str1,1)
			elif "o1" in smiles:
				print "case 5"
				ion_smiles = smiles
				#replace_str1 = "[" + ion + "](o1)"
				replace_str1 = "o1([" + ion + "])"
				ion_smiles = ion_smiles.replace("o1", replace_str1,1)
			elif "O" in smiles:
 				print "case 6"
 				ion_smiles = "" + smiles + ""
 				#replace_str1 = "[" + ion + "]O"
 				replace_str1 = "O([" + ion + "])"
 				ion_smiles = ion_smiles.replace("O", replace_str1,1)
			elif "C#N" in smiles:
				print "case 7"
				ion_smiles = smiles
				replace_str1 = "C#N([" + ion + "])"
				ion_smiles = ion_smiles.replace("C#N", replace_str1,1)
			elif "n1" in smiles:
				print "case 8"
				ion_smiles = smiles
				replace_str1 = "n1([" + ion + "])"
				ion_smiles = ion_smiles.replace("n1", replace_str1,1)
			elif "n2" in smiles:
				print "case 9"
				ion_smiles = smiles
				replace_str1 = "n2([" + ion + "])"
				ion_smiles = ion_smiles.replace("n2", replace_str1,1)
			elif "S" in smiles:
				print "case 10"
				ion_smiles = smiles
				replace_str1 = "S([" + ion + "])"
				ion_smiles = ion_smiles.replace("S", replace_str1,1)
			elif "s1" in smiles:
				print "case 11"
				ion_smiles = smiles
				replace_str1 = "s1([" + ion + "])"
				ion_smiles = ion_smiles.replace("s1", replace_str1,1)
			else:
				print "case 12"
				ion_smiles = "(" + smiles + ")"
				ion_smiles = ion_smiles + "[" + ion + "]"
	
			print "ion_smiles", ion_smiles

			tmpxyz = workdir + "/tmp.xyz"
			babel_str = 'obabel ' + '-:"' + ion_smiles + '" -oxyz -O ' + tmpxyz + ' -h --gen3D'
			os.system("%s >& /dev/null" % babel_str)

			ion_solv  = read(tmpxyz)
		#
		# end making initial structure for metal-solvent complex
		#

		#
		# set charge
		#
		if ion in ["Li", "Na", "K", "Rb", "Cs"]:
			ion_charge = 1
		elif ion in ["Mg", "Be", "Ca", "Sr", "Ba"]:
			ion_charge = 2
		else:
			print "only alkali and alkaline earth metal allowed" ; quit()

	else:
		print "found in database"
		ion_charge = db_ion.get(num=num).calculator_parameters["charge"]
		ion_smiles = db_ion.get(num=num).smiles
		delete_num_from_json(num, solv_ion_jsonfile)

	if "gau" in calculator:
 		label = workdir + "/g16_low_" + ion
 		ion_solv.calc = Gaussian(label=label, method=xc, basis=basis, scf=scf, charge=ion_charge, mult=1, gfinput="gfinput", grid=grid, nprocshared=nprocs, mem=mem)
		traj = ion + "_low" + str(num).zfill(4) + ".traj"
 		FIRE(ion_solv, trajectory=traj).run(fmax=fmax,steps=steps)
		## avoid using Gaussian because it stops at non-convergence
 		#ion_solv.calc = Gaussian(method=xc, basis=basis, label=label, opt=opt, force=None,
		# 						 charge=ion_charge, mult=1, nprocshared=12, ioplist=ioplist)
		#ion_solv.get_potential_energy()
	elif "nw" in calculator:
 		label = workdir + "/nwchem_low_" + ion
	 	ion_solv.calc = NWChem(label=label, xc=xc, basis=basis, charge=ion_charge, mult=1, 
							   iterations=200, mulliken=True, memory=memory) # cation
                        
 		traj = ion + "_low_" + str(num).zfill(4) + ".traj"
 		FIRE(ion_solv).run(fmax=fmax)

 	db_ion.write(ion_solv, smiles=ion_smiles, name=name, num=num, 
				molecular_weight=wgt, density=dens,
 				boiling_point=bp, melting_point=mp, flushing_point=fp, pubchemCID=pubchem,
 			  	level="low" )

shutil.rmtree(workdir)

