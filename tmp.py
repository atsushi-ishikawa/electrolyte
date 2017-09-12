from ase import Atoms
from mynwchem import NWChem

h2 = Atoms("H2",[(0,0,0),(0,0,0.7)])
calc = NWChem(xc="b3lyp",mulliken=True)
h2.set_calculator(calc)

h2.get_potential_energy()

c1 = calc.get_homo_center()
c2 = calc.get_lumo_center()

mul = calc.results["mul_charge"]

print c1, c2
print mul
