from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

calc = GPAW(mode=PW(500),
            xc='PBE',
            h=0.10,
            occupations=FermiDirac(width=0.001),
            kpts={'size': (8, 8, 8), 'gamma': True},
            parallel={'band': 1, 'domain': 1},
            txt='gs_Bi2Te3.txt')

a = 4.3835
c = 30.487
mu = 0.400
nu = 0.212
cell = [[-a / 2, -3**0.5 / 6 * a, c / 3],
        [a / 2, -3**0.5 / 6 * a, c / 3],
        [0.0, 3**0.5 / 3 * a, c / 3]]
pos = [[mu, mu, mu],
       [-mu, -mu, -mu],
       [0.0, 0.0, 0.0],
       [nu, nu, nu],
       [-nu, -nu, -nu]]
a = Atoms('Bi2Te3', cell=cell, scaled_positions=pos, pbc=True)
a.set_calculator(calc)
a.get_potential_energy()

calc.write('gs_Bi2Te3.gpw')
