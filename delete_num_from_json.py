from ase.db import connect
import sys

argvs = sys.argv
db = connect(argvs[1])
id = db.get(num=int(argvs[2])).id
db.delete([id])

