"""
Compares Source/Matter/DiracFieldCartoon.hpp with the Python reference (dirac_ref.cartoon_rhs / cartoon_sources) at
random points. Build and run:
    g++ -O2 -std=c++14 -DCH_SPACEDIM=2 -Istub -I../../../Source/Matter -I../../../Source/utils test_dirac_cpp.cpp -o test_dirac_cpp
    python3 test_dirac_cpp.py ./test_dirac_cpp
"""
import os, sys, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'derivations'))
import numpy as np
from dirac_ref import cartoon_rhs, cartoon_sources

exe = sys.argv[1] if len(sys.argv) > 1 else './test_dirac_cpp'
rng = np.random.default_rng(11)
lines, refs = [], []
for trial in range(200):
    sym = lambda s: (lambda m: (m + m.T) / 2)(rng.normal(size=(2, 2)) * s)
    h = np.eye(2) + sym(0.2)
    v = dict(chi=rng.uniform(0.3, 1.2), h=h, hww=rng.uniform(0.8, 1.2), K=rng.normal() * 0.3, A=sym(0.2),
             Aww=rng.normal() * 0.2, lapse=rng.uniform(0.3, 1.0), shift=rng.normal(size=2) * 0.2)
    v['dchi'] = list(rng.normal(size=2) * 0.3)
    v['dh'] = [sym(0.3) for k in range(2)]
    v['dhww'] = list(rng.normal(size=2) * 0.3)
    v['dlapse'] = list(rng.normal(size=2) * 0.3)
    v['dshift'] = rng.normal(size=(2, 2)) * 0.2            # [k, i] = d_k beta^i
    v['psi'] = rng.normal(size=4) + 1j * rng.normal(size=4)
    v['dpsi'] = [rng.normal(size=4) + 1j * rng.normal(size=4) for k in range(2)]
    y, mu = rng.uniform(0.05, 3.), rng.uniform(0.5, 1.5)
    rhs, aux = cartoon_rhs(mu, v, y)
    adv = sum(aux['beta'][k] * aux['dpsi'][k] for k in range(2))
    rho, Si, Sij, Sww, n = cartoon_sources(mu, v, y)
    refs.append(np.concatenate([np.ravel(np.column_stack([(rhs - adv).real, (rhs - adv).imag])),
                                [rho, Si[0], Si[1], Sij[0, 0], Sij[0, 1], Sij[1, 1], Sww, n]]))
    vals = [v['chi'], h[0, 0], h[0, 1], h[1, 1], v['hww'], v['K'], v['A'][0, 0], v['A'][0, 1], v['A'][1, 1], v['Aww'],
            v['lapse'], v['shift'][0], v['shift'][1]]
    for k in range(2):
        vals += [v['dchi'][k], v['dh'][k][0, 0], v['dh'][k][0, 1], v['dh'][k][1, 1], v['dhww'][k], v['dlapse'][k],
                 v['dshift'][k, 0], v['dshift'][k, 1]]
    vals += list(v['psi'].real) + list(v['psi'].imag)
    vals += [v['dpsi'][k][c].real for c in range(4) for k in range(2)]
    vals += [v['dpsi'][k][c].imag for c in range(4) for k in range(2)]
    vals += [y, mu]
    lines.append(' '.join('%.17e' % x for x in vals))
out = subprocess.run([exe], input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout
cpp = np.array([[float(x) for x in l.split()] for l in out.strip().split('\n')])
ref = np.array(refs)
err = np.max(np.abs(cpp - ref) / (1 + np.abs(ref)))
print('C++ vs Python reference at %d random points: max relative error %.1e  %s' % (len(ref), err, 'OK' if err < 1e-11 else 'FAILED'))
