import os
import glob
import commands
from ase.io import read
from ase.db import connect

dirlist = glob.glob("/fs01/electrolyte/[0-9][0-9]*")

inlist = []

for i in dirlist:
	tmp = glob.glob(i + "/*.log")
	if len(tmp) != 0:
		inlist.append(tmp[0])

con = connect("tmp.json")

for i in inlist
	babel_str1 = "obabel " + " -ig03 " + inlist[i] + " -ocan "
	tmp = commands.getoutput(babel_str1)
	smiles = tmp.split("\t")[0]

	babel_str2 = "obabel " + " -ig03 " + inlist[i] + " -oxyz -O tmp.xyz"
	os.system(babel_str2)
	
	mol = read("tmp.xyz")
	con.write(mol,smiles=smiles)

