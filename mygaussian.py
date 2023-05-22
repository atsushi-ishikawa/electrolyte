def prepare_basisfile(atoms, basisfile=None, basisname=None, ecplist=None):
	import os

	basis_homedir="/home/a_ishi/ase/basis_and_ecp_for_gaussian"
	basisname = basisname.lower()

	elems = atoms.get_chemical_symbols()
	elems = set(elems)

	f_tot = open(basisfile,"w")

	for elem in elems:
		file = elem + "_" + basisname + ".bas"
		file = os.path.join(basis_homedir, file)

		if os.path.isfile(file):
			f_elem = open(file, "r")
		else:
			print("error opening basis %s file for %s" % (basisname,elem))
			exit()

		lines = f_elem.readlines()
		f_tot.writelines(lines)
		f_elem.close()

	f_tot.write("\n")

	if ecplist is None:
		pass
	else:
		for elem in ecplist:
			file = elem + "_" + ecplist[elem] + ".ecp"
			file = os.path.join(basis_homedir, file)

			if os.path.isfile(file):
				f_elem = open(file, "r")
			else:
				print("error opening ECP %s file for %s" % (basisname,elem))
				exit()

			lines = f_elem.readlines()
			f_tot.writelines(lines)
			f_elem.close()

	f_tot.write("\n")
	f_tot.close()

