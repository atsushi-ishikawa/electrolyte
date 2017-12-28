## procedure to do electrolyte calc

### 1. Calculation for solvent

* When Atoms in ASE class is ready, do calculation by `solv.py`. Results will be stored in JSON file specified in this script.
* When Atoms is not ready, do small calculation by `prepare_atoms.py` and use resultant JSON file for solvent calculation.
* Original JSON file should be
```python
{
    "0": {
        "boiling_point": 241.6,
        "density": 1.2,
        "flushing_point": 135.0,
        "melting_point": -48.8,
        "pubchemCID": 7924,
        "smiles": "CC1COC(=O)O1",
        "molecular_weight": 102.089,
        "name": "PROPYLENE CARBONATE",
        "num": 1
    }
    ...
}
```
* `solv.py 1` where 1 is the solvent number (`num` key in JSON file)
* If `solv.py` is finished, new JSON file is formed.


#### Performing calculation
One by one
```bash
qsub run.sh
```
Multiple
```bash
python submit.py 1 19
```

### 2. Make guess for ion coordinated solvents
* Edit `ion.py` accordingly
* `ion.py Li Na 1` where 1 is the solvent number
* New JSON file is generated. Can be the same with input JSON file.

### 3. Final calculation for ion coordinated solvents
* Edit `coordinated.py` accordingly
* `coordinated.py Li Na 1` where 1 is the solvent number
* Final JSON file is generated

### 4. Converting JSON to CSV
* Copy JSON file from server if necessary
```bash {cmd=True}
scp whisky:/home/a_ishi/ase/electrolyte/ishi3_final.json ./
```
* Convert JSON to CSV for Excel use
```python {cmd=True}
from tools import json_to_csv
json_to_csv("ishi3_final.json","tmp.csv")
```

* Drop unnecessary keys and wrong calculations
```python {cmd=True}
import pandas
df = pandas.read_csv('tmp.csv')

# delete unnecessary keys (to regression)
drop_col = ['calculator', 'charge', 'iterations', 'memory',
  'mulliken', 'mult', 'xc', 'cell', 'ctime',
  'R_ion_O.Li', 'R_ion_O.Na', 'R_ion_O.K', 'R_ion_O.Rb','R_ion_O.Cs',
  'e_homo_ion.Li', 'e_homo_ion.Na',   'e_homo_ion.K',
  'e_homo_ion.Rb', 'e_homo_ion.Cs',
  'e_lumo_ion.Li', 'e_lumo_ion.Na',   'e_lumo_ion.K',
  'e_lumo_ion.Rb', 'e_lumo_ion.Cs',
  'mul_ion_ion.Li', 'mul_ion_ion.Na', 'mul_ion_ion.K',
  'mul_ion_ion.Rb', 'mul_ion_ion.Cs',
  'low_ion_ion.Li', 'low_ion_ion.Na', 'low_ion_ion.K',
  'low_ion_ion.Rb', 'low_ion_ion.Cs',
  'mul_O_ion.Li',   'mul_O_ion.Na',   'mul_O_ion.K',
  'mul_O_ion.Rb',   'mul_O_ion.Cs',
  'low_O_ion.Li',   'low_O_ion.Na',   'low_O_ion.K',
  'low_O_ion.Rb',   'low_O_ion.Cs',
  'user', 'unique_id', 'positions', 'pbc', 'numbers', 'mtime', 'magmom',
  'dipole', 'energy', 'forces'
  ]
df = df.drop(drop_col, axis=1)

# delete erroneous calculations
df = df[ df['Ecoord.Li'] > -3.0]

df.to_csv('tmp.csv')
```

* Plot data by seaborn
```python {cmd=True}
import pandas
import matplotlib.pyplot as plt
import seaborn as sns

df = pandas.read_csv("tmp.csv")
# sns.jointplot("low_O_solv","Ecoord.Li",data=df)

sns.set()
sns.set(font_scale=0.8)
sns.pairplot(df, vars=[
  'Ecoord.Li','mul_O_solv','low_O_solv','e_homo','e_lumo',
  'iso_polarizability',
  'total_dipole', 'ionization_potential', 'electron_affinity'
  ], plot_kws={"s": 15}, diag_kind="hist", size=1.4)
plt.show()
```
