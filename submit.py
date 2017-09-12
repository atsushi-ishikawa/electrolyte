import os,sys

argvs = sys.argv
num_from = int(argvs[1])
num_to   = int(argvs[2])

for i in range(num_from, num_to+1):
	command = "qsub run.sh " + str(i)
	os.system(command)

