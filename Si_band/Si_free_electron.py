from ase.build import bulk
from ase.calculators.test import FreeElectrons

a = bulk('Si', 'diamond', 5.43)
a.calc = FreeElectrons(nvalence=1,
                       kpts={'path': 'WLGXWK', 'npoints': 200})
a.get_potential_energy()
bs = a.calc.band_structure()
bs.plot(emax=30.0, filename='Si_free.png')