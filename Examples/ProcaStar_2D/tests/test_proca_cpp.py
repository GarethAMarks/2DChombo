"""
Compares Source/Matter/ProcaFieldCartoon.hpp with the Python reference (proca_cartoon_ref) at random points.
Build and run:
    g++ -O2 -std=c++14 -DCH_SPACEDIM=2 -Istub -I../../../Source/Matter -I../../../Source/utils test_proca_cpp.cpp -o test_proca_cpp
    python3 test_proca_cpp.py ./test_proca_cpp
"""
import os, sys, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'derivations'))
import numpy as np
from proca_cartoon_ref import proca_cartoon_rhs_and_sources

exe = sys.argv[1] if len(sys.argv) > 1 else './test_proca_cpp'
rng = np.random.default_rng(23)
cv = lambda *s: rng.normal(size=s) + 1j * rng.normal(size=s)
lines, refs = [], []
for trial in range(200):
    sym = lambda s: (lambda m: (m + m.T) / 2)(rng.normal(size=(2, 2)) * s)
    chi, h, hww, Kt = rng.uniform(0.3, 1.2), np.eye(2) + sym(0.2), rng.uniform(0.8, 1.2), rng.normal() * 0.3
    lapse, shift = rng.uniform(0.3, 1.0), rng.normal(size=2) * 0.2
    dchi, dh, dhww = rng.normal(size=2) * 0.3, [sym(0.3) for k in range(2)], rng.normal(size=2) * 0.2
    dlapse, dshift = rng.normal(size=2) * 0.3, rng.normal(size=(2, 2)) * 0.2       # dshift[k, i] = d_k beta^i
    X, E, P, Z = cv(2), cv(2), cv(1)[0], cv(1)[0]
    dX, dE, dP, dZ = cv(2, 2), cv(2, 2), cv(2), cv(2)                             # [k, i] = d_k
    ddX = cv(2, 2, 2); ddX = 0.5 * (ddX + np.transpose(ddX, (1, 0, 2)))           # [j, k, i] = d_j d_k X_i
    y, mu, kappa = rng.uniform(0.05, 3.), rng.uniform(0.5, 1.5), rng.uniform(0., 2.)
    g2, gww = h / chi, hww / chi
    dg2 = np.array([dh[k] / chi - h * dchi[k] / chi**2 for k in range(2)])
    dgww = dhww / chi - hww * dchi / chi**2
    # any trace-free A_ij: only the trace K enters the Proca equations
    K2, Kww = h * Kt / 3 / chi, hww * Kt / 3 / chi
    rhs, src = proca_cartoon_rhs_and_sources(mu, kappa, y, lapse, dlapse, shift, dshift, g2, dg2, gww, dgww, K2, Kww,
                                             dict(X=X, E=E, Xphi=P, Z=Z), dict(X=dX, E=dE, Xphi=dP, Z=dZ), ddX)
    rX, rE = rhs['X'] - shift @ dX, rhs['E'] - shift @ dE
    rP, rZ = rhs['Xphi'] - shift @ dP, rhs['Z'] - shift @ dZ
    refs.append(np.concatenate([rX.real, rX.imag, rE.real, rE.imag, [rP.real, rP.imag, rZ.real, rZ.imag, src['rho']],
                                src['Si'], [src['Sij'][0, 0], src['Sij'][0, 1], src['Sij'][1, 1], src['Sww'], src['n'],
                                            src['Gauss'].real, src['Gauss'].imag]]))
    vals = [chi, h[0, 0], h[0, 1], h[1, 1], hww, Kt, lapse, shift[0], shift[1]]
    for k in range(2):
        vals += [dchi[k], dh[k][0, 0], dh[k][0, 1], dh[k][1, 1], dhww[k], dlapse[k], dshift[k, 0], dshift[k, 1]]
    f12 = lambda X, E, P, Z: list(X.real) + list(X.imag) + list(E.real) + list(E.imag) + [P.real, P.imag, Z.real, Z.imag]
    vals += f12(X, E, P, Z)
    for k in range(2):
        vals += f12(dX[k], dE[k], dP[k], dZ[k])
    for j in range(2):
        for k in range(2):
            vals += list(ddX[j, k].real) + list(ddX[j, k].imag)
    vals += [y, mu, kappa]
    lines.append(' '.join('%.17e' % x for x in vals))
out = subprocess.run([exe], input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout
cpp = np.array([[float(x) for x in l.split()] for l in out.strip().split('\n')])
ref = np.array(refs)
err = np.max(np.abs(cpp - ref) / (1 + np.abs(ref)))
print('C++ vs Python reference at %d random points: max relative error %.1e  %s' % (len(ref), err, 'OK' if err < 1e-11 else 'FAILED'))
