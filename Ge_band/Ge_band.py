"""Band structure

Calculate the band structure of Ge along high symmetry directions
Brillouin zone
"""

import numpy as np
from ase.build import bulk
from ase.dft.kpoints import ibz_points, get_bandpath
from gpaw import GPAW, PW, FermiDirac

# Perform standard ground state calculation (with plane wave basis)
ge = bulk('Ge', 'diamond', 5.6575)
calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 8, 8),
            random=True,  # random guess (needed if many onoccupied bands required)
            occupations=FermiDirac(0.01),
            txt='Ge_gs.txt')
ge.calc = calc
ge.get_potential_energy()
calc.write('Ge_gs.gpw')
ef = calc.get_fermi_level()

nbands = 8

# Restart from ground state and fix potential:
calc = GPAW('Ge_gs.gpw',
            nbands=16,
            fixdensity=True,
            symmetry='off',
            convergence={'bands': nbands})

# Use ase.dft module for obtaining k-points along high symmetry directions
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
L = points['L']
kpts, x, X = get_bandpath([W, L, G, X, W, K], calc.atoms.cell, npoints=60)
calc.set(kpts=kpts)
calc.get_potential_energy()
e_kn = np.array([calc.get_eigenvalues(k) for k in range(len(kpts))])

# Plot the band structure

bs = calc.band_structure()
bs.plot(filename='Ge_band.png', show=True, emax=15.0)

import matplotlib.pyplot as plt

e_kn -= ef
emin = e_kn.min() - 1.0
emax = e_kn[:, nbands].max() + 1.0

plt.figure(figsize=(5, 6))
for n in range(nbands):
    plt.plot(x, e_kn[:, n])
for p in X:
    plt.plot([p, p], [emin, emax], 'k-')
plt.plot([0, X[-1]], [0, 0], 'k-')
plt.xticks(X, ['$%s$' % n for n in ['W', 'L', r'\Gamma', 'X', 'W', 'K']])
plt.axis(xmin=0, xmax=X[-1], ymin=emin, ymax=emax)
plt.xlabel('k-vector')
plt.ylabel('E - E$_F$ (eV)')
plt.title('Bandstructure of Germanium')
plt.savefig('Ge_band_colored.png')
