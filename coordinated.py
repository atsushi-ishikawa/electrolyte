from ase import Atoms
from mynwchem import NWChem
from ase.db import connect
from ase.optimize import BFGS
from ase.io import read, write
from tools import ABcoord
import os
import sys
import numpy as np
from numpy import linalg
import mynwchem

ions = []
argvs = sys.argv

nions = len(argvs) - 2 

IP_and_EA = True

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
basis    = "3-21G"
fmax     =  0.10
memory   = "total 8 gb"
response = "1 0.0"
# ---

# - label -
label_solv = "calc" + str(num).zfill(4) + "/nwchem_solv"
# ---------

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
charge = db_solv.get(num=num).calculator_parameters["charge"]
mult   = db_solv.get(num=num).calculator_parameters["mult"]

print " --- calculating num = ",num, "(", name,")"
print "solvent"
#  --- neutral ---
#
# geometry optimization
#
calc = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge, 
              mult=mult, iterations=200, mulliken=True, memory=memory)
solv.set_calculator(calc)
traj = "solv_" + str(num).zfill(4) + ".traj"
BFGS(solv, trajectory=traj).run(fmax=fmax, steps=maxsteps)
#
# single point calculation for polarizability
#
calc = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge, mult=mult,
	      iterations=200, mulliken=True, memory=memory, response=response)
solv.set_calculator(calc)
E_solv = solv.get_potential_energy()
#
# --- analyze results ---
#
mul_charge = calc.results["mul_charge"] 
low_charge = calc.results["low_charge"] 

e_homo = calc.get_homo_energy()
e_lumo = calc.get_lumo_energy()

polar     = calc.results["polarizability"]
iso_pol   = polar[0]
aniso_pol = polar[1]

dipole       = calc.get_dipole_moment() ; dipole = np.array(dipole)
total_dipole = linalg.norm(dipole)
#
# Get mulliken charge of solvent's O atom
# which is "going to" coordinate to ion.
# Coordinating O is determined from highest MO 
# with large coef on O.
#
(mo_cent, element, occ) = calc.get_mo_center()

O_cent = calc.get_atomcentered_occupied_mo("O")
highest_O_cent_mo = O_cent[-1]
O_cent_mo = mo_cent[highest_O_cent_mo]

mul_O_solv = mul_charge[O_cent_mo-1]
low_O_solv = low_charge[O_cent_mo-1]
#
#  --- IP and EA calculation ---
#
if IP_and_EA:
	#
	#  cation
	#
	calc_c = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge+1, 
			mult=mult+1, iterations=200, mulliken=True, memory=memory)
	solv.set_calculator(calc_c)
	traj = "solv_cat" + str(num).zfill(4) + ".traj"
	BFGS(solv, trajectory=traj).run(fmax=fmax, steps=maxsteps)
	E_cat = solv.get_potential_energy()
	#
	#  anion
	#
	calc_a = NWChem(label=label_solv, xc=xc, basis=basis, charge=charge-1,
			mult=mult+1, iterations=200, mulliken=True, memory=memory)
	solv.set_calculator(calc_a)
	traj = "solv_ani" + str(num).zfill(4) + ".traj"
	BFGS(solv, trajectory=traj).run(fmax=fmax, steps=maxsteps)
	E_ani = solv.get_potential_energy()
	#
	IP = E_cat - E_solv
	EA = E_ani - E_solv
else:
	IP = 0.0; EA = 0.0

# -----------------------------------
# --------- ion coordinated ---------
# -----------------------------------
e_homo_ion  = {} ; e_lumo_ion  = {}
mul_ion_ion = {} ; mul_O_ion   = {}
low_ion_ion = {} ; low_O_ion   = {}
Ecoord  = {}     ; R_ion_O = {}

for ion in ions:
	ion_jsonfile = solv_jsonfile.split(".")[0] + "_" + ion + ".json"
	db_ion       = connect(ion_jsonfile)

	label_ion = "calc" + ion + str(num).zfill(4) + "/nwchem"

	ion_solv = db_ion.get_atoms(num=num)
	charge   = db_ion.get(num=num).calculator_parameters["charge"]
	mult     = db_ion.get(num=num).calculator_parameters["mult"]

	print ion, " coordinated"
	calc = NWChem(label=label_ion, xc=xc, basis=basis, charge=charge, 
                      mult=mult, iterations=200, 
   	              mulliken=True, memory=memory) # cation
	ion_solv.set_calculator(calc)

	traj = ion + str(num).zfill(4) + ".traj"
	BFGS(ion_solv, trajectory=traj).run(fmax=fmax, steps=maxsteps)

	E_ion_solv = ion_solv.get_potential_energy()

	#
	# Mulliken and Lowdin charge
	#
	mul_charge_ion = calc.results["mul_charge"]
	low_charge_ion = calc.results["low_charge"]

	e_homo_ion[ion] = calc.get_homo_energy()
	e_lumo_ion[ion] = calc.get_lumo_energy()

	R_ion_O[ion], O_idx = ABcoord(ion_solv, ion, "O")

	lst = ion_solv.get_chemical_symbols()
	ion_idx = lst.index(ion)

	mul_ion_ion[ion] = mul_charge_ion[ion_idx] # ion
	mul_O_ion[ion]   = mul_charge_ion[O_idx]   # coordinating O of O-Li
	low_ion_ion[ion] = low_charge_ion[ion_idx] # ion
	low_O_ion[ion]   = low_charge_ion[O_idx]   # coordinating O of O-Li

	#
	# get energy of ion 
	#
	if ion in ["Li", "Na", "K", "Rb", "Cs"]:
		ion_charge = 1
	elif ion == "Mg":
		ion_charge = 2
	else:
		print "only alkaline metal and Mg is allowed" ; quit()

	label_atom = "calc" + ion + "/nwchem"
	ion_atom = Atoms(ion, positions=[(0,0,0)])
	ion_atom.calc = NWChem(label=label_atom, xc=xc, basis=basis, 
			       charge=ion_charge, mulliken=True)

	E_ion = ion_atom.get_potential_energy()

	Ecoord[ion] = E_ion_solv - ( E_ion + E_solv )

# end ion loop

db_out.write(solv,num=num, smiles=smiles, name=name,
			  	molecular_weight=wgt, density=dens,
			  	boiling_point=bp, melting_point=mp, flushing_point=fp,
				data={	"e_homo"       : e_homo,
					"e_lumo"       : e_lumo,
					"e_homo_ion"   : e_homo_ion,
					"e_lumo_ion"   : e_lumo_ion,
					"mul_O_solv"   : mul_O_solv,
					"mul_ion_ion"  : mul_ion_ion,
					"mul_O_ion"    : mul_O_ion,
					"low_O_solv"   : low_O_solv,
					"low_ion_ion"  : low_ion_ion,
					"low_O_ion"    : low_O_ion,
					"R_ion_O"      : R_ion_O,
					"total_dipole" : total_dipole,
					"iso_polarizability"   : iso_pol,
					"aniso_polarizability" : aniso_pol,
					"ionization_potential" : IP,
					"electron_affinity"    : EA,
					"Ecoord"      : Ecoord }
			)
