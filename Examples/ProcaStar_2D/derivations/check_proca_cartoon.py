"""
Checks the cartoon reduction of the Proca equations (proca_cartoon_ref.py) against the full 3D equations
(proca3d_ref.py, itself checked against the covariant 4D theory and the 1D code by check_proca3d.py).
Run: python3 check_proca_cartoon.py

Arbitrary twist-free axisymmetric data (metric, lapse, shift, K_ij, complex X_i, E^i, X_phi, Z built from functions of
(x, rho) in cylindrical form) are evaluated in 3D with finite differences in all three directions, and the 3D right-hand
side and sources at (x, y, 0) are compared with the cartoon ones computed from x and y derivatives only. Also checks
that the azimuthal components of the 3D right-hand side (X_z, E^z) and of S_i, S_xz, S_yz vanish on the plane, so the
reduced system is closed.
"""
import numpy as np
from proca3d_ref import proca_rhs_and_sources
from proca_cartoon_ref import proca_cartoon_rhs_and_sources

rng = np.random.default_rng(11)
c = rng.normal(size=80)
mu, kappa = 0.7, 0.4
h = 4e-3


def D_(f, X, k, hh=h):
    e = np.zeros(len(X)); e[k] = hh
    return sum(w * (f(X + j * e) - f(X - j * e)) for j, w in [(1, 3 / 4), (2, -3 / 20), (3, 1 / 60)]) / hh


def sm(n, x, r, amp=0.2):   # smooth function of (x, rho^2)
    return amp * np.sin(c[n] + c[n + 1] * x + c[n + 2] * r * r)


def fields(p):
    """twist-free axisymmetric fields at the 3D point p: cylindrical (x, rho) components -> Cartesian"""
    x, y, z = p
    r = np.hypot(y, z); ey = np.array([0., y / r, z / r]); ex = np.array([1., 0, 0]); ep = np.array([0., -z / r, y / r])
    a, b, cc, d = 1 + sm(0, x, r), r * sm(3, x, r, 0.1), 1 + sm(6, x, r), 1 + sm(9, x, r)
    gamma = a * np.outer(ex, ex) + b * (np.outer(ex, ey) + np.outer(ey, ex)) + cc * np.outer(ey, ey) + d * np.outer(ep, ep)
    a, b, cc, d = sm(12, x, r), r * sm(15, x, r), sm(18, x, r), sm(21, x, r)
    K = a * np.outer(ex, ex) + b * (np.outer(ex, ey) + np.outer(ey, ex)) + cc * np.outer(ey, ey) + d * np.outer(ep, ep)
    alpha = 1 + sm(24, x, r)
    beta = sm(27, x, r) * ex + r * sm(30, x, r) * ey
    cx = lambda n: (0.5 + sm(n, x, r)) * np.exp(1j * (c[n + 3] + 0.5 * x + c[n + 4] * r * r))
    X = cx(33) * ex + r * cx(38) * ey
    E = cx(43) * ex + r * cx(48) * ey
    return dict(alpha=alpha, beta=beta, gamma=gamma, K=K, X=X, E=E, Xphi=cx(53), Z=0.3 * cx(58))


ok = True
for p0 in [np.array([0.3, 0.7, 0.]), np.array([-0.5, 0.25, 0.]), np.array([0.1, 1.6, 0.])]:
    f0 = fields(p0)
    d3 = lambda key, k: D_(lambda q: fields(q)[key], p0, k)
    dalpha = np.array([d3('alpha', k) for k in range(3)]); dbeta = np.array([d3('beta', k) for k in range(3)])
    dgamma = np.array([d3('gamma', k) for k in range(3)])
    dV = {key: np.array([d3(key, k) for k in range(3)]) for key in ['X', 'E', 'Xphi', 'Z']}
    ddX = np.array([[D_(lambda q: D_(lambda w: fields(w)['X'], q, l), p0, k) for l in range(3)] for k in range(3)])
    V = {key: f0[key] for key in ['X', 'E', 'Xphi', 'Z']}
    rhs3, src3 = proca_rhs_and_sources(mu, kappa, f0['alpha'], dalpha, f0['beta'], dbeta, f0['gamma'], dgamma, f0['K'],
                                       V, dV, ddX)
    # cartoon: in-plane components and x, y derivatives only
    s2 = slice(0, 2)
    V2 = dict(X=f0['X'][s2], E=f0['E'][s2], Xphi=f0['Xphi'], Z=f0['Z'])
    dV2 = dict(X=dV['X'][s2, s2], E=dV['E'][s2, s2], Xphi=dV['Xphi'][s2], Z=dV['Z'][s2])
    rhs2, src2 = proca_cartoon_rhs_and_sources(
        mu, kappa, p0[1], f0['alpha'], dalpha[s2], f0['beta'][s2], dbeta[s2, s2], f0['gamma'][s2, s2],
        dgamma[s2, s2, s2], f0['gamma'][2, 2], dgamma[s2, 2, 2], f0['K'][s2, s2], f0['K'][2, 2], V2, dV2,
        ddX[s2, s2, s2])
    err = max(np.max(np.abs(rhs2['X'] - rhs3['X'][s2])), np.max(np.abs(rhs2['E'] - rhs3['E'][s2])),
              abs(rhs2['Xphi'] - rhs3['Xphi']), abs(rhs2['Z'] - rhs3['Z']))
    errs = max(abs(src2['rho'] - src3['rho']), np.max(np.abs(src2['Si'] - src3['Si'][s2])),
               np.max(np.abs(src2['Sij'] - src3['Sij'][s2, s2])), abs(src2['Sww'] - src3['Sij'][2, 2]),
               abs(src2['Gauss'] - src3['Gauss']), abs(src2['n'] - src3['n']))
    zero = max(abs(rhs3['X'][2]), abs(rhs3['E'][2]), abs(src3['Si'][2]), abs(src3['Sij'][0, 2]), abs(src3['Sij'][1, 2]))
    for name, e in [('RHS (X_x, X_y, E^x, E^y, X_phi, Z)', err), ('sources (rho, S_i, S_ij, S_ww), Gauss, n', errs),
                    ('azimuthal components of the 3D RHS and sources', zero)]:
        good = e < 1e-8
        ok &= good
        print('  at (x, y) = (%5.2f, %4.2f): %-50s %.1e  %s' % (p0[0], p0[1], name, e, 'OK' if good else 'FAILED'))
print('\nALL CHECKS PASSED' if ok else '\nSOME CHECKS FAILED')
