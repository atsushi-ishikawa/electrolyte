from tools import json_to_csv
import os
import sys

argvs = sys.argv

json_to_csv(str(argvs[1]),"tmp.csv")

os.system("open tmp.csv")

