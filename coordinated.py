from ase import Atoms
from mynwchem import NWChem
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from ase.optimize import FIRE
from ase.io import read, write
from tools import ABcoord
import os, sys, shutil
import numpy as np
from numpy import linalg
import mynwchem

ions = []
argvs = sys.argv

nions = len(argvs) - 2 

IP_and_EA = False
polarizability = False

calculator = "gaussian"
calculator = calculator.lower()

population = "nbo"

# get ion list
for i in range(0, nions):
	ions.append(str(argvs[i+1]))

num   = int(argvs[nions + 1])

#solv_jsonfile = "electrolye_2017Aug.json"
#out_jsonfile  = "tmpout.json"
solv_jsonfile = "ishi3_new.json"
out_jsonfile  = "ishi3_final.json"

# ---
xc       = "B3LYP"
#basis    = "cc-pvdz"
basis    = "svp"
fmax     =  0.10
memory   = "total 8 gb"
response = "1 0.0"
# ---

# workdir
work = "/work/a_ishi/"
if not os.path.isdir(work):
	os.makedirs(work)

workdir = work + "calc" + str(num).zfill(4)

# label
label_solv = workdir + "/solv"

db_solv = connect(solv_jsonfile)
db_out  = connect(out_jsonfile)

maxsteps = 200 # maximum number of geometry optimization steps

# ---------------------------
# --------- solvent ---------
# ---------------------------
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

# get charge and multiplicity from old file
try:
	charge = db_solv.get(num=num).calculator_parameters["charge"]
	mult   = db_solv.get(num=num).calculator_parameters["mult"]
except:
	charge = 0
	mult   = 1

print " --- calculating num = ",num, "(", name,")"
print "solvent"
#  --- neutral ---
#
# geometry optimization
#
if "gau" in calculator:
 	calc = Gaussian(label=label_solv, method=xc, basis=basis, charge=charge, mult=mult)
elif "nw" in calcullator:
	calc = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge, mult=mult, iterations=200, mulliken=True, memory=memory)

solv.set_calculator(calc)
traj = "solv_" + str(num).zfill(4) + ".traj"
FIRE(solv, trajectory=traj).run(fmax=fmax, steps=maxsteps)
#
# single point calculation
#
if "gau" in calculator:
 	calc = Gaussian(label=label_solv, method=xc, basis=basis, charge=charge, mult=mult, force=None, pop="full,nbo")
elif "nw" in calculator:
	if polarizability:
		calc = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge, mult=mult,
					  iterations=200, mulliken=True, memory=memory, response=response)
	else:
		calc = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge, mult=mult,
					  iterations=200, mulliken=True, memory=memory)

solv.set_calculator(calc)
E_solv = solv.get_potential_energy()
#
# --- analyze results ---
#
if "gau" in calculator:
	atomic_charge = solv.get_nbo_charge()
elif "nw" in calculator:
	if population=="mulliken":
		atomic_charge = calc.results["mul_charge"] 
	elif population=="lowdin":
		atomic_charge = calc.results["low_charge"] 
	else:
		print "choose population"
		exit()

if "gau" in calculator:
	mo_energy = solv.get_mo_energy()
	e_homo = mo_energy[0][-1]
	e_lumo = mo_energy[1][0]
elif "nw" in calculator:
	e_homo = calc.get_homo_energy()
	e_lumo = calc.get_lumo_energy()

if polarizability:
	polar = calc.results["polarizability"]
	iso_pol   = polar[0]
	aniso_pol = polar[1]

dipole       = calc.get_dipole_moment() ; dipole = np.array(dipole)
total_dipole = linalg.norm(dipole)
#
# Get atomic charge of solvent's O atom which is "going to" coordinate to ion.
# Coordinating O is determined from highest MO with large coef on O.
#
if "gau" in calculator:
	mocoef, basis_to_atom = solv.get_mo_coeff()
	nocc    = len(mo_energy[0])
	c_homo  = np.array( mocoef[nocc-1] )
	maxind  = np.argmax(map(abs, c_homo))
	O_coord = basis_to_atom[maxind] # atom index where HOMO coefficient is maximum --> coordinating O atom
	O_solv_charge = atomic_charge[O_coord]
elif "nw" in calculator:
	(mo_cent, element, occ) = calc.get_mo_center()
	O_coord = calc.get_atomcentered_occupied_mo("O")
	highest_O_coord_mo = O_coord[-1]
	O_coord_mo = mo_cent[highest_O_coord_mo]
	O_solv_charge = atomic_charge[O_coord_mo - 1]
#
#  --- IP and EA calculation ---
#
if IP_and_EA:
	#
	#  cation
	#
	print "calculating cation"
	solv_c = solv.copy()
	if "gau" in calculator:
		calc_c = Gaussian(label=label_solv, method=xc, basis=basis, charge=chg, mult=mlt, population="full,nbo")
	else:
		calc_c = NWChem(label=label_solv, xc=xc, basis=basis, charge=chg, mult=mlt, iterations=200, mulliken=True, memory=memory)

	chg = solv_c.get_initial_charges()
	mlt = solv_c.get_initial_magnetic_moments()

	chg[0] += 1 ; mlt[0] += 1

	solv_c.set_initial_charges(charges=chg)
	solv_c.set_initial_magnetic_moments(magmoms=mlt)
	solv_c.set_calculator(calc_c)

	traj = "solv_cat" + str(num).zfill(4) + ".traj"
	FIRE(solv_c, trajectory=traj).run(fmax=fmax, steps=maxsteps)
	E_cat = solv_c.get_potential_energy()
	#
	#  anion
	#
	print "calculating anion"
	solv_a = solv.copy()
	if "gau" in calculator:
		calc_a = Gaussian(label=label_solv, method=xc, basis=basis, charge=chg, mult=mlt, population="full,nbo")
	else:
		calc_a = NWChem(label=label_solv, xc=xc, basis=basis, charge=chg, mult=mlt, iterations=200, mulliken=True, memory=memory)

	chg = solv_a.get_initial_charges()
	mlt = solv_a.get_initial_magnetic_moments()

	chg[0] -= 1 ; mlt[0] += 1

	solv_a.set_initial_charges(charges=chg)
	solv_a.set_initial_magnetic_moments(magmoms=mlt)
	solv_a.set_calculator(calc_a)

	traj = "solv_ani" + str(num).zfill(4) + ".traj"
	FIRE(solv_a, trajectory=traj).run(fmax=fmax, steps=maxsteps)
	E_ani = solv_a.get_potential_energy()
	#
	IP = E_cat - E_solv
	EA = E_ani - E_solv
else:
	IP = 0.0; EA = 0.0

# -----------------------------------
# --------- ion coordinated ---------
# -----------------------------------
#  descriptors for coordinates system
#    coordination energy (Ecoord) and ion-O distance (R_ion_O)
#
Ecoord  = {}
R_ion_O = {}

for ion in ions:
	ion_jsonfile = solv_jsonfile.split(".")[0] + "_" + ion + ".json"
	db_ion       = connect(ion_jsonfile)

	label_ion = workdir + "/coord_" + ion

	ion_solv = db_ion.get_atoms(num=num)
	try:
		charge = db_ion.get(num=num).calculator_parameters["charge"]
		mult   = db_ion.get(num=num).calculator_parameters["mult"]
	except:
		charge = 0
		mult   = 1

	print ion, " coordinated"
	if "gau" in calculator:
		calc = Gaussian(label=label_ion, method=xc, basis=basis, charge=charge, mult=mult, population="full,nbo") # cat
	else:
		calc = NWChem(label=label_ion, xc=xc, basis=basis, charge=charge, 
					  mult=mult, iterations=200, mulliken=True, memory=memory) # cat
	ion_solv.set_calculator(calc)

	traj = ion + str(num).zfill(4) + ".traj"
	FIRE(ion_solv, trajectory=traj).run(fmax=fmax, steps=maxsteps)

	E_ion_solv = ion_solv.get_potential_energy()
	R_ion_O[ion], O_idx = ABcoord(ion_solv, ion, "O")

	lst = ion_solv.get_chemical_symbols()
	ion_idx = lst.index(ion)
	#
	# get energy of ion 
	#
	if ion in ["Li", "Na", "K", "Rb", "Cs"]:
		ion_charge = 1
	elif ion == "Mg":
		ion_charge = 2
	else:
		print "only alkaline metal and Mg is allowed"
		exit()

	ion_atom = Atoms(ion, positions=[(0,0,0)])
	if "gau" in calculator:
		label_atom = "calc" + ion + "/g16"
		ion_atom.calc = Gaussian(label=label_atom, method=xc, basis=basis, 
								 charge=ion_charge, population="full,nbo")
	elif "nw" in calculator:
		label_atom = "calc" + ion + "/nwchem"
		ion_atom.calc = NWChem(label=label_atom, xc=xc, basis=basis, 
							   charge=ion_charge, mulliken=True)

	E_ion = ion_atom.get_potential_energy()

	Ecoord[ion] = E_ion_solv - ( E_ion + E_solv )

	print "Ecoord for",ion," = ",Ecoord[ion]
	print E_ion_solv, E_ion, E_solv

# end ion loop

db_out.write(solv,num=num, smiles=smiles, name=name,
			  	molecular_weight=wgt, density=dens,
			  	boiling_point=bp, melting_point=mp, flushing_point=fp,
				data={	"e_homo"        : e_homo,
						"e_lumo"        : e_lumo,
						"O_solv_charge" : O_solv_charge,
						"R_ion_O"       : R_ion_O,
						"total_dipole"  : total_dipole,
					#	"iso_polarizability"   : iso_pol,
					#	"aniso_polarizability" : aniso_pol,
					#	"ionization_potential" : IP,
					#	"electron_affinity"    : EA,
						"Ecoord"        : Ecoord }
			)

#shutil.rmtree(workdir)

