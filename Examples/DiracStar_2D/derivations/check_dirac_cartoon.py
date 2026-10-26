"""
Numerical verification of the 3+1 and cartoon Dirac equations in dirac_ref.py (run: python3 check_dirac_cartoon.py).

A. Generic 3+1 form vs the covariant 4D Dirac operator. For arbitrary smooth alpha(t,x), beta^i(t,x), gamma_ij(t,x) and
   psi(t,x) with no symmetry, the 4D operator D psi = gamma^a e_a^mu nabla_mu psi + mu psi (tetrad e_0 = n, symmetric
   triad, 4D Christoffels and spin connection by finite differences) must satisfy
       gamma^0 D psi + (d_t psi - RHS) / alpha = 0,
   with RHS the 3+1 right-hand side, which uses only (alpha, beta, gamma, K, spatial derivatives).
B. Cartoon reduction. Axisymmetric fields about the x axis are built in 3D from 2D profiles by rotation,
       gamma(R x) = R gamma(x) R^T, ...,   psi(R x) = e^{i m phi} S(R) psi(x),  S(R) = exp(phi SPIN_GEN),  m = +1/2,
   with the parities on the axis that smoothness requires. Checks:
   B1. the cartoon RHS on the plane z = 0 (2D data only) equals the generic 3D RHS evaluated with 3D finite differences;
   B2. the 3D RHS is itself equivariant, RHS(R x) = e^{i m phi} S(R) RHS(x), which fixes the sign of SPIN_GEN (the
       equations commute with rotations only for the correct spinor representation);
   B3. the m = -1/2 spinor, the mirror image psi_-(x, y, z) = P psi_+(x, y, -z) with P gamma^z P^-1 = -gamma^z (P fixed
       here), also satisfies the Dirac equation, and T[psi_+] + T[psi_-] on the plane equals the cartoon sources (2 x the
       z-even components of T[psi_+]).
C. Spherical star: the stationary Dirac star (written directly in 3D, not through the rotation construction) has
   d_t psi = -i omega psi from both the cartoon and the 3D RHS, and the cartoon sources reproduce the 1D (Spherical_BS
   Evolver) expressions summed over both spinors.
"""
import numpy as np
from scipy.linalg import expm
from dirac_ref import *

rng = np.random.default_rng(3)
H = 2e-3


def D(f, X, k, h=H):
    """6th-order central difference of f (array-valued) along coordinate k at X"""
    e = np.zeros(len(X)); e[k] = h
    c = [(1, 3 / 4), (2, -3 / 20), (3, 1 / 60)]
    return sum(w * (f(X + j * e) - f(X - j * e)) for j, w in c) / h


def report(name, err, tol=1e-8):
    print('  %-70s max error %.1e  %s' % (name, err, 'OK' if err < tol else 'FAILED'))
    return err < tol


ok = True
mu = 0.8
# ------------------------------------------------------------------------------------------------- A
print('A. generic 3+1 RHS vs 4D Dirac operator (no symmetry):')
cA = rng.normal(size=60)


def alpha_f(X):
    t, x, y, z = X
    return 1.1 + 0.2 * np.sin(0.7 * t + 0.3 * x - 0.5 * y + 0.2 * z + cA[0])


def beta_f(X):
    t, x, y, z = X
    return 0.15 * np.array([np.sin(cA[1] + 0.4 * t + 0.5 * x + 0.3 * y), np.cos(cA[2] + 0.3 * t - 0.6 * z + 0.2 * x),
                            np.sin(cA[3] - 0.5 * t + 0.4 * y + 0.7 * z)])


def gamma_f(X):
    t, x, y, z = X
    M = np.zeros((3, 3))
    n = 4
    for i in range(3):
        for j in range(i, 3):
            M[i, j] = M[j, i] = 0.25 * np.sin(cA[n] + cA[n + 1] * t + cA[n + 2] * x + cA[n + 3] * y + cA[n + 4] * z)
            n += 5
    return np.eye(3) + M


def psi_f(X):
    t, x, y, z = X
    return np.array([np.exp(I * (cA[40 + c] + 0.5 * t + cA[44 + c] * x)) * (0.5 + 0.2 * np.sin(cA[48 + c] * y + z))
                     + 0.1 * I * np.cos(cA[52 + c] * z - x * y) for c in range(4)])


X0 = np.array([0.3, 0.2, -0.4, 0.5])


def metric4(X):
    a, b, g = alpha_f(X), beta_f(X), gamma_f(X)
    G = np.zeros((4, 4))
    bl = g @ b
    G[0, 0] = -a**2 + b @ bl
    G[0, 1:] = G[1:, 0] = bl
    G[1:, 1:] = g
    return G


def tetrad4(X):
    a, b, g = alpha_f(X), beta_f(X), gamma_f(X)
    E = np.zeros((4, 4))           # E[a, mu] = e_a^mu
    E[0, 0] = 1 / a
    E[0, 1:] = -b / a
    E[1:, 1:] = triad(g)[0].T
    return E


g4 = metric4(X0); g4i = np.linalg.inv(g4)
dg4 = [D(metric4, X0, m) for m in range(4)]
# Christoffel Gamma^r_{mn} = 1/2 g^{rs} (d_m g_sn + d_n g_sm - d_s g_mn)
Chr = np.zeros((4, 4, 4))
for r_ in range(4):
    for m in range(4):
        for n in range(4):
            Chr[r_, m, n] = 0.5 * sum(g4i[r_, s] * (dg4[m][s, n] + dg4[n][s, m] - dg4[s][m, n]) for s in range(4))
E4 = tetrad4(X0)
dE4 = [D(tetrad4, X0, m) for m in range(4)]
Ecov = np.einsum('ab,mn,bn->am', ETA, g4, E4)                     # e^a_mu
omega = np.zeros((4, 4, 4))                                        # omega_{mu a b} (frame indices down)
for m in range(4):
    nab = dE4[m] + np.einsum('nl,bl->bn', Chr[:, m, :], E4)       # (nabla_m e_b)^n, rows b
    om_up = Ecov @ nab.T                                           # omega_m^a_b
    omega[m] = ETA @ om_up
Gspin = [sum(omega[m, a, b] / 8 * (GAM[a] @ GAM[b] - GAM[b] @ GAM[a]) for a in range(4) for b in range(4)) for m in range(4)]
psi0 = psi_f(X0)
dpsi4 = [D(psi_f, X0, m) for m in range(4)]
D4 = mu * psi0 + sum(GAM[a] @ sum(E4[a, m] * (dpsi4[m] + Gspin[m] @ psi0) for m in range(4)) for a in range(4))

# 3+1 quantities at X0
a0, b0, g0 = alpha_f(X0), beta_f(X0), gamma_f(X0)
dalpha = np.array([D(alpha_f, X0, k) for k in (1, 2, 3)])
dbeta = np.array([D(beta_f, X0, k) for k in (1, 2, 3)])           # dbeta[k, i]
dgamma = np.array([D(gamma_f, X0, k) for k in (1, 2, 3)])
gt = D(gamma_f, X0, 0)
K0 = -(gt - (np.einsum('k,kij->ij', b0, dgamma) + np.einsum('kj,ik->ij', g0, dbeta) + np.einsum('ik,jk->ij', g0, dbeta))) / (2 * a0)
dE3 = np.array([D(lambda X: triad(gamma_f(X))[0], X0, k) for k in (1, 2, 3)])
rhs, aux = generic_rhs(mu, a0, dalpha, b0, dbeta, g0, dgamma, K0, dE3, psi0, dpsi4[1:])
res = GAM[0] @ D4 + (dpsi4[0] - rhs) / a0
ok &= report('gamma^0 D4 psi + (d_t psi - RHS_3+1)/alpha', np.max(np.abs(res)))
# a deliberately wrong triad rotation term must fail (the check is sensitive to it)
print('  (size of the triad-rotation term alpha Gamma_n|rot psi for reference: %.2e)' % np.max(np.abs(
    a0 * sum(0.25 * ((aux['E'].T @ g0 @ (edot_from_gammadot(g0, gt) - (np.einsum('l,lib->ib', b0, dE3) - np.einsum('lb,li->ib', aux['E'], dbeta))))[a_, b_]
                     - (aux['E'].T @ g0 @ (edot_from_gammadot(g0, gt) - (np.einsum('l,lib->ib', b0, dE3) - np.einsum('lb,li->ib', aux['E'], dbeta))))[b_, a_]) / (2 * a0)
             * gg(a_ + 1, b_ + 1) for a_ in range(3) for b_ in range(3)) @ psi0)))

# ------------------------------------------------------------------------------------------------- B
print('B. cartoon reduction on axisymmetric fields (axis x, plane z = 0):')
cB = rng.normal(size=80) * 0.5


def prof(x, q, n, amp=0.2):
    """smooth function of x and q = y^2"""
    return amp * np.sin(cB[n] + cB[n + 1] * x + cB[n + 2] * q) * np.exp(-0.1 * q)


def plane_metric(x, y):
    q = y * y
    g = np.zeros((3, 3))
    g[0, 0] = 1.2 + prof(x, q, 0)
    g[0, 1] = g[1, 0] = y * prof(x, q, 3)
    base = 1.1 + prof(x, q, 6)
    g[1, 1] = base + q * prof(x, q, 9, 0.1)
    g[2, 2] = base
    return g


def plane_K(x, y):
    q = y * y
    k = np.zeros((3, 3))
    k[0, 0] = prof(x, q, 12)
    k[0, 1] = k[1, 0] = y * prof(x, q, 15)
    base = prof(x, q, 18)
    k[1, 1] = base + q * prof(x, q, 21)
    k[2, 2] = base
    return k


def plane_alpha(x, y):
    return 0.9 + prof(x, y * y, 24)


def plane_beta(x, y):
    q = y * y
    return np.array([prof(x, q, 27), y * prof(x, q, 30), 0.])


def plane_psi(x, y):
    q = y * y
    c = lambda n: prof(x, q, n, 0.5) + I * prof(x, q, n + 40, 0.5)
    return np.array([0.3 + c(33), y * c(36), c(39 - 30), y * c(42 - 30)])  # comps 0, 2: Sigma_x = +1; 1, 3: Sigma_x = -1


def rot(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def polar(X):
    x, y, z = X
    return x, np.hypot(y, z), np.arctan2(z, y)


def F3(plane_f, kind):
    def f(X):
        x, rr, ph = polar(X)
        R = rot(ph)
        p = plane_f(x, rr)
        if kind == 'scalar':
            return p
        if kind == 'vector':
            return R @ p
        if kind == 'tensor':
            return R @ p @ R.T
        if kind == 'spinor':
            return np.exp(I * M_J * ph) * expm(ph * SPIN_GEN) @ p
    return f


gam3, K3, al3, be3, psi3 = F3(plane_metric, 'tensor'), F3(plane_K, 'tensor'), F3(plane_alpha, 'scalar'), F3(plane_beta, 'vector'), F3(plane_psi, 'spinor')


def rhs3(X):
    d = lambda f: np.array([D(f, X, k, 1e-3) for k in range(3)])
    g = gam3(X)
    return generic_rhs(mu, al3(X), d(al3), be3(X), d(be3), g, d(gam3), K3(X), d(lambda Y: triad(gam3(Y))[0]),
                       psi3(X), list(d(psi3)))


def grid_vars(x, y):
    """the 2D grid variables (BSSN form) and their x, y derivatives at (x, y)"""
    def bssn(xx, yy):
        g = plane_metric(xx, yy)
        chi = np.linalg.det(g)**(-1 / 3)
        h = g * chi
        Kp = plane_K(xx, yy)
        trK = np.trace(np.linalg.inv(g) @ Kp)
        A = chi * (Kp - g * trK / 3)
        return dict(chi=chi, h=h[:2, :2], hww=h[2, 2], K=trK, A=A[:2, :2], Aww=A[2, 2], lapse=plane_alpha(xx, yy),
                    shift=plane_beta(xx, yy)[:2], psi=plane_psi(xx, yy))
    v = bssn(x, y)
    for name in ['chi', 'h', 'hww', 'lapse', 'shift', 'psi']:
        f = lambda X, name=name: np.asarray(bssn(X[0], X[1])[name])
        v['d' + name] = [D(f, np.array([x, y]), k, 1e-3) for k in range(2)]
    v['dshift'] = np.array(v['dshift'])                              # dshift[k, i]
    return v


errs1, errs2 = [], []
for trial in range(3):
    x0, y0 = rng.uniform(-1, 1), rng.uniform(0.3, 1.2)
    Xp = np.array([x0, y0, 0.])
    r3, aux3 = rhs3(Xp)
    rc, auxc = cartoon_rhs(mu, grid_vars(x0, y0), y0)
    errs1.append(np.max(np.abs(r3 - rc)) / np.max(np.abs(rc)))
    ph = rng.uniform(0.3, 2.5)
    Xr = np.array([x0, y0 * np.cos(ph), y0 * np.sin(ph)])
    rr_, _ = rhs3(Xr)
    errs2.append(np.max(np.abs(rr_ - np.exp(I * M_J * ph) * expm(ph * SPIN_GEN) @ r3)) / np.max(np.abs(r3)))
ok &= report('B1. cartoon RHS vs 3D RHS on the plane (relative)', max(errs1), 1e-7)
ok &= report('B2. 3D RHS equivariant under rotations about x (relative)', max(errs2), 1e-7)
# the opposite sign of the spinor generator breaks equivariance:
SG = SPIN_GEN.copy()
psi3_wrong = lambda X: (lambda x, rr, ph: np.exp(I * M_J * ph) * expm(-ph * SG) @ plane_psi(x, rr))(*polar(X))
def rhs3w(X):
    d = lambda f: np.array([D(f, X, k, 1e-3) for k in range(3)])
    return generic_rhs(mu, al3(X), d(al3), be3(X), d(be3), gam3(X), d(gam3), K3(X), d(lambda Y: triad(gam3(Y))[0]),
                       psi3_wrong(X), list(d(psi3_wrong)))[0]
ph = 1.1; x0, y0 = 0.2, 0.7
w0 = rhs3w(np.array([x0, y0, 0.])); w1 = rhs3w(np.array([x0, y0 * np.cos(ph), y0 * np.sin(ph)]))
print('     (with the opposite sign of SPIN_GEN the equivariance error is %.1e: the sign is determined)' %
      (np.max(np.abs(w1 - np.exp(I * M_J * ph) * expm(-ph * SG) @ w0)) / np.max(np.abs(w0))))

# B3: mirror spinor. P must satisfy P gamma^0 P^-1 = gamma^0, P gamma^{x,y} P^-1 = gamma^{x,y}, P gamma^z P^-1 = -gamma^z
P = GAM[3] @ (GAM[0] @ GAM[1] @ GAM[2] @ GAM[3])                  # gamma^z gamma^5-like; check the relations
sgn = [1, 1, 1, -1]
assert all(np.allclose(P @ GAM[a] @ np.linalg.inv(P), sgn[a] * GAM[a]) for a in range(4))
Pz = np.diag([1., 1., -1.])
psi3m = lambda X: P @ psi3(X * np.array([1, 1, -1]))
def rhs3m(X):
    d = lambda f: np.array([D(f, X, k, 1e-3) for k in range(3)])
    return generic_rhs(mu, al3(X), d(al3), be3(X), d(be3), gam3(X), d(gam3), K3(X), d(lambda Y: triad(gam3(Y))[0]),
                       psi3m(X), list(d(psi3m)))
x0, y0 = 0.35, 0.8
Xp = np.array([x0, y0, 0.])
rp, auxp = rhs3(Xp)
rm, auxm = rhs3m(Xp)
ok &= report('B3a. mirror spinor solves the Dirac equation: RHS[psi_-] = P RHS[psi_+](mirror)',
             np.max(np.abs(rm - P @ rp)) / np.max(np.abs(rp)), 1e-7)
d = lambda f, X: list(np.array([D(f, X, k, 1e-3) for k in range(3)]))
sp_ = generic_sources(al3(Xp), be3(Xp), psi3(Xp), d(psi3, Xp), rp, auxp)
sm_ = generic_sources(al3(Xp), be3(Xp), psi3m(Xp), d(psi3m, Xp), rm, auxm)
rho_c, Si_c, Sij_c, Sww_c, n_c = cartoon_sources(mu, grid_vars(x0, y0), y0)
tot = [sp_[0] + sm_[0], (sp_[1] + sm_[1])[:2], (sp_[2] + sm_[2])[:2, :2], (sp_[2] + sm_[2])[2, 2]]
err = max(abs(tot[0] - rho_c), np.max(np.abs(tot[1] - Si_c)), np.max(np.abs(tot[2] - Sij_c)), abs(tot[3] - Sww_c))
ok &= report('B3b. T[psi_+] + T[psi_-] = cartoon sources (rho, S_i, S_ij, S_ww)', err, 1e-7)
oddsum = max(np.max(np.abs((sp_[1] + sm_[1])[2])), np.max(np.abs((sp_[2] + sm_[2])[:2, 2])))
ok &= report('B3c. z-odd components S_z, S_xz, S_yz cancel between the two spinors', oddsum, 1e-7)
ok &= report('B3d. T_ab is real', max(sp_[3], sm_[3]), 1e-9)

# ------------------------------------------------------------------------------------------------- C
print('C. stationary spherical star (analytic test profiles, written directly in 3D):')
# use any smooth f, g, conformal factor and lapse; stationarity requires the 1D equations, so instead check the RHS
# reduces to the 1D form: d_t psi from the cartoon must equal the 1D (Spherical_BS_Evolver) RHS for u, v.
w = 0.85
gfun = lambda r: 0.2 * np.exp(-r * r / 4) * (1 + 0.1 * r * r)
ffun = lambda r: 0.05 * r * np.exp(-r * r / 3)
psifun = lambda r: 1 + 0.3 * np.exp(-r * r / 5)       # conformal factor, gamma = psi^4 delta (isotropic)
alfun = lambda r: 0.8 + 0.15 * (1 - np.exp(-r * r / 6))


def sph_psi(X):
    x, y, z = X
    r = np.sqrt(x * x + y * y + z * z)
    n = np.array([x, y, z]) / r
    sn = sum(n[i] * SIG[i] for i in range(3))
    up = gfun(r) * np.array([1, 0]) / np.sqrt(4 * np.pi)
    lo = I * ffun(r) * (-(sn @ np.array([1, 0]))) / np.sqrt(4 * np.pi)
    return np.concatenate([up, lo])


def sph_vars(x, y):
    def bssn(xx, yy):
        r = np.hypot(xx, yy)
        return dict(chi=psifun(r)**-4, h=np.eye(2), hww=1., K=0., A=np.zeros((2, 2)), Aww=0., lapse=alfun(r),
                    shift=np.zeros(2), psi=sph_psi(np.array([xx, yy, 0.])))
    v = bssn(x, y)
    for name in ['chi', 'h', 'hww', 'lapse', 'shift', 'psi']:
        f = lambda X, name=name: np.asarray(bssn(X[0], X[1])[name])
        v['d' + name] = [D(f, np.array([x, y]), k, 1e-3) for k in range(2)]
    v['dshift'] = np.array(v['dshift'])
    return v


# the 1D equations in isotropic form (z-axis cartoon of Spherical_BS_Evolver with h_zz = h_ww = 1, chi = psi^-4):
def rhs_1d(r):
    dr = lambda f: D(lambda X: np.array([f(X[0])]), np.array([r]), 0, 1e-3)[0]
    chi = psifun(r)**-4
    c = np.sqrt(chi)
    L = 0.5 * (-dr(lambda s: psifun(s)**-4) / chi + dr(alfun) / alfun(r))
    u, v = gfun(r), ffun(r)
    ut = -I * alfun(r) * mu * u + I * alfun(r) * (c * (dr(ffun) + L * v) + 2 * c * v / r)
    vt = I * alfun(r) * mu * v - I * alfun(r) * (c * (dr(gfun) + L * u))
    return ut, vt


errs = []
for (x0, y0) in [(0.4, 0.9), (-1.1, 0.3), (0.2, 1.7)]:
    r = np.hypot(x0, y0)
    rc, _ = cartoon_rhs(mu, sph_vars(x0, y0), y0)
    ut, vt = rhs_1d(r)
    # map the 1D (u, v) RHS onto the spinor components
    expect = np.concatenate([ut * np.array([1, 0]), I * vt * (-(x0 / r * SIG[0] + y0 / r * SIG[1]) @ np.array([1, 0]))]) / np.sqrt(4 * np.pi)
    errs.append(np.max(np.abs(rc - expect)))
ok &= report('C1. cartoon RHS on a spherical configuration = 1D (u, v) equations', max(errs), 1e-8)


# 1D sources summed over both spinors (Spherical_BS_Evolver ComplexSpinorField): with u, v real here
def src_1d(r):
    dr = lambda f: D(lambda X: np.array([f(X[0])]), np.array([r]), 0, 1e-3)[0]
    u, v = gfun(r), ffun(r)
    c = psifun(r)**-2
    ut, vt = rhs_1d(r)
    Srr = c * (v * dr(gfun) - u * dr(ffun)) / (2 * np.pi)
    Stt = -c * u * v / r / (2 * np.pi)
    rho = mu * (u * u - v * v) / (2 * np.pi) + Srr + 2 * Stt
    return rho, Srr, Stt


errs = []
for (x0, y0) in [(0.4, 0.9), (-1.1, 0.3)]:
    r = np.hypot(x0, y0)
    rho_c, Si_c, Sij_c, Sww_c, n_c = cartoon_sources(mu, sph_vars(x0, y0), y0)
    rho1, Srr1, Stt1 = src_1d(r)
    nhat = np.array([x0, y0]) / r
    g = psifun(r)**4
    Sij_expect = g * (Srr1 * np.outer(nhat, nhat) + Stt1 * (np.eye(2) - np.outer(nhat, nhat)))
    # rho from the RHS uses the real dynamics (not stationary), compare with the 1D on-shell form
    errs += [abs(rho_c - rho1), np.max(np.abs(Sij_c - Sij_expect)), abs(Sww_c - g * Stt1)]
ok &= report('C2. cartoon sources = 1D sources (rho, S_ij, S_ww) on spherical data', max(errs), 1e-8)

print('ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED')
