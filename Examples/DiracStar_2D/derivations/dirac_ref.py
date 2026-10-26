"""
Reference (Python) implementation of the 3+1 Dirac equations used by DiracField.hpp, in two forms:

  generic_rhs / generic_sources : any 3-metric, lapse, shift and extrinsic curvature in 3D (no symmetry assumed).
  cartoon_rhs / cartoon_sources : the same equations on the cartoon plane z = 0 of an axisymmetric spacetime
                                  (symmetry axis x, cylindrical radius y), using only the 2D grid data. This is a
                                  line-by-line mirror of the C++ in Source/Matter/DiracField.impl.hpp.

Conventions (as in Spherical_BS_Evolver's DiracStar, signature (-,+,+,+)):
  * (gamma^a nabla_a + mu) psi = 0, {gamma^a, gamma^b} = 2 eta^{ab}, gamma^a = -i gamma_D^a with gamma_D the standard
    Dirac representation built from the Pauli matrices permuted so that the spin quantisation axis is the symmetry axis
    x:  (sigma_x, sigma_y, sigma_z) -> (sigma_3, sigma_1, sigma_2).
  * psibar = psi^dag beta,  beta = i gamma^0 = gamma_D^0.
  * Spinor connection Gamma_mu = (1/8) omega_{mu ab} [gamma^a, gamma^b],  omega_mu^a_b = e^a_nu nabla_mu e_b^nu.
  * T_ab = (1/4)(X_ab + X_ba),  X_ab = psibar gamma_a nabla_b psi - (nabla_a psi)^dag beta gamma_b psi.
  * Tetrad: e_0 = n (time gauge) and the symmetric spatial triad e_a^i = (gamma^{-1/2})^{ia} (Cartesian components),
    an algebraic function of the 3-metric that is covariant under rotations, hence smooth on the axis.
  * K_ij = -nabla_i n_j (so d_t gamma_ij = -2 alpha K_ij + L_beta gamma_ij), as in GRChombo.

3+1 form (derived in check_dirac_cartoon.py):
  d_t psi = beta^k d_k psi - alpha Gamma_n psi + alpha gamma^0 gamma^b e_b^k (d_k psi + Gamma_k psi) + alpha mu gamma^0 psi,
  Gamma_k = (1/2) omega_{k0b} gamma^0 gamma^b + (1/4) omega_{kab} gamma^a gamma^b,   omega_{k0b} = K_kj e_b^j,
  Gamma_n = (1/2) n.omega_{0b} gamma^0 gamma^b + (1/4) n.omega_{ab} gamma^a gamma^b,  n.omega_{0b} = -e_b^j d_j alpha / alpha,
  n.omega_{ab} = (1/(2 alpha)) [ e_a . (d_t e_b - L_beta e_b) - (a <-> b) ]  (rotation of the triad; d_t e_b follows from
                 d_t gamma_ij through the symmetric square root),
  omega_{kab} = gamma_ij e_a^i (d_k e_b^j + Gamma^j_kl e_b^l)  (spatial spin connection).
"""
import numpy as np
from scipy.linalg import sqrtm, solve_sylvester

I = 1j
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [s3, s1, s2]          # sigma'_x, sigma'_y, sigma'_z: spin quantised along the symmetry axis x
Z2, I2 = np.zeros((2, 2), dtype=complex), np.eye(2, dtype=complex)


def block(a, b, c, d):
    return np.block([[a, b], [c, d]])


gD = [block(I2, Z2, Z2, -I2)] + [block(Z2, s, -s, Z2) for s in SIG]
GAM = [-I * g for g in gD]                                    # gamma^a (upper frame index)
ETA = np.diag([-1., 1., 1., 1.])
GAM_LOW = [ETA[a, a] * GAM[a] for a in range(4)]              # gamma_a
BETA = I * GAM[0]                                             # = gamma_D^0
SIGX = block(s3, Z2, Z2, s3)                                  # Sigma_x: spin along the axis
M_J = 0.5                                                     # the evolved spinor has J_x = +1/2


def gg(a, b):
    return GAM[a] @ GAM[b]


# --------------------------------------------------------------------------------------------------------------------
# generic 3D
# --------------------------------------------------------------------------------------------------------------------
def triad(gamma):
    """e_a^i as the columns of E = gamma^{-1/2}; co-triad theta = gamma^{1/2} (theta^a_i = theta[a, i])"""
    E = np.real(sqrtm(np.linalg.inv(gamma)))
    return E, np.linalg.inv(E)


def christoffel(gamma, dgamma):
    gi = np.linalg.inv(gamma)
    G = np.zeros((3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                G[i, j, k] = 0.5 * sum(gi[i, l] * (dgamma[j, l, k] + dgamma[k, l, j] - dgamma[l, j, k]) for l in range(3))
    return G


def gamma_dot(alpha, beta, dbeta, gamma, dgamma, K):
    """d_t gamma_ij = -2 alpha K_ij + beta^k d_k gamma_ij + gamma_kj d_i beta^k + gamma_ik d_j beta^k; dbeta[k, i] = d_k beta^i"""
    gd = -2 * alpha * K + np.einsum('k,kij->ij', beta, dgamma)
    gd += np.einsum('kj,ik->ij', gamma, dbeta) + np.einsum('ik,jk->ij', gamma, dbeta)
    return gd


def connections(alpha, dalpha, beta, dbeta, gamma, dgamma, K, E, dE, Edot):
    """returns Gamma_k (k = 0..2, spatial coordinate directions) and Gamma_n (4x4 matrices)"""
    G3 = christoffel(gamma, dgamma)
    # spatial omega_{k ab}, a, b spatial frame indices 1..3 -> 0..2 here
    om = np.zeros((3, 3, 3))
    for k in range(3):
        nab_e = dE[k] + np.einsum('jl,lb->jb', G3[:, k, :], E)          # (nabla_k e_b)^j as columns
        om[k] = E.T @ gamma @ nab_e
    om0 = np.einsum('kj,jb->kb', K, E)                                   # omega_{k0b} = K_kj e_b^j
    # time parts
    an0 = -(E.T @ dalpha) / alpha                                        # n.omega_{0b}
    Lbe = np.einsum('l,lib->ib', beta, dE) - np.einsum('lb,li->ib', E, dbeta)   # (L_beta e_b)^i
    W = E.T @ gamma @ (Edot - Lbe)
    anab = (W - W.T) / (2 * alpha)                                       # n.omega_{ab}
    Gk = []
    for k in range(3):
        M = np.zeros((4, 4), dtype=complex)
        for b in range(3):
            M += 0.5 * om0[k, b] * gg(0, b + 1)
            for a in range(3):
                M += 0.25 * om[k, a, b] * gg(a + 1, b + 1)
        Gk.append(M)
    Gn = np.zeros((4, 4), dtype=complex)
    for b in range(3):
        Gn += 0.5 * an0[b] * gg(0, b + 1)
        for a in range(3):
            Gn += 0.25 * anab[a, b] * gg(a + 1, b + 1)
    return Gk, Gn


def edot_from_gammadot(gamma, gd):
    """d_t of E = sqrt(gamma^{-1}): solve E X + X E = d_t(gamma^{-1})"""
    gi = np.linalg.inv(gamma)
    E = np.real(sqrtm(gi))
    return solve_sylvester(E, E, -gi @ gd @ gi)


def generic_rhs(mu, alpha, dalpha, beta, dbeta, gamma, dgamma, K, dE, psi, dpsi):
    E, th = triad(gamma)
    gd = gamma_dot(alpha, beta, dbeta, gamma, dgamma, K)
    Edot = edot_from_gammadot(gamma, gd)
    Gk, Gn = connections(alpha, dalpha, beta, dbeta, gamma, dgamma, K, E, dE, Edot)
    rhs = sum(beta[k] * dpsi[k] for k in range(3)) - alpha * Gn @ psi + alpha * mu * GAM[0] @ psi
    for b in range(3):
        nb = sum(E[k, b] * (dpsi[k] + Gk[k] @ psi) for k in range(3))
        rhs += alpha * gg(0, b + 1) @ nb
    return rhs, dict(E=E, th=th, Gk=Gk, Gn=Gn)


def frame_derivs(alpha, beta, psi, dpsi, psi_t, E, Gk, Gn):
    """nabla_a psi for a = 0 (normal) and the spatial frame directions"""
    nab = [(psi_t - sum(beta[k] * dpsi[k] for k in range(3))) / alpha + Gn @ psi]
    for b in range(3):
        nab.append(sum(E[k, b] * (dpsi[k] + Gk[k] @ psi) for k in range(3)))
    return nab


def stress_frame(psi, nab):
    psibar = psi.conj() @ BETA
    X = np.zeros((4, 4), dtype=complex)
    for a in range(4):
        for b in range(4):
            X[a, b] = psibar @ GAM_LOW[a] @ nab[b] - nab[a].conj() @ BETA @ GAM_LOW[b] @ psi
    T = 0.25 * (X + X.T)
    return T


def generic_sources(alpha, beta, psi, dpsi, psi_t, aux):
    """rho, S_i (= -gamma_i^a n^b T_ab), S_ij for ONE spinor; also returns the imaginary part of T as a check"""
    T = stress_frame(psi, frame_derivs(alpha, beta, psi, dpsi, psi_t, aux['E'], aux['Gk'], aux['Gn']))
    th = aux['th']
    rho = T[0, 0]
    S_i = -th.T @ T[1:, 0]
    S_ij = th.T @ T[1:, 1:] @ th
    return rho.real, S_i.real, S_ij.real, np.max(np.abs(T.imag))


# --------------------------------------------------------------------------------------------------------------------
# cartoon (mirror of the C++)
# --------------------------------------------------------------------------------------------------------------------
GEN = np.array([[0., 0., 0.], [0., 0., -1.], [0., 1., 0.]])     # d R_phi / d phi at phi = 0 (rotation about x)
SPIN_GEN = -0.5 * GAM[2] @ GAM[3]                               # spinor generator of that rotation, = -(i/2) Sigma_x
                                                                # (sign fixed by check B2 in check_dirac_cartoon.py)


def sqrt2_and_deriv(M, dMs):
    """S = sqrt(M) for SPD 2x2 M, S = (M + s I)/t, s = sqrt(det M), t = sqrt(tr M + 2 s); and dS for each dM"""
    s = np.sqrt(np.linalg.det(M))
    t = np.sqrt(np.trace(M) + 2 * s)
    S = (M + s * np.eye(2)) / t
    Minv = np.linalg.inv(M)
    out = []
    for dM in dMs:
        ds = 0.5 * s * np.trace(Minv @ dM)
        dt = (np.trace(dM) + 2 * ds) / (2 * t)
        out.append((dM + ds * np.eye(2)) / t - S * dt / t)
    return S, out


def cartoon_fields(v, y):
    """3x3 metric quantities on the plane from the 2D grid variables (dict v holding values and x,y derivatives)"""
    chi, h, hww = v['chi'], v['h'], v['hww']
    gamma = np.zeros((3, 3))
    gamma[:2, :2] = h / chi
    gamma[2, 2] = hww / chi
    dgamma = np.zeros((3, 3, 3))
    for k in range(2):
        dgamma[k, :2, :2] = v['dh'][k] / chi - h * v['dchi'][k] / chi**2
        dgamma[k, 2, 2] = v['dhww'][k] / chi - hww * v['dchi'][k] / chi**2
    dgamma[2] = (GEN @ gamma - gamma @ GEN) / y
    K = np.zeros((3, 3))
    K[:2, :2] = (v['A'] + h * v['K'] / 3) / chi
    K[2, 2] = (v['Aww'] + hww * v['K'] / 3) / chi
    beta = np.array([v['shift'][0], v['shift'][1], 0.])
    dbeta = np.zeros((3, 3))
    dbeta[:2, :2] = v['dshift']                                   # dbeta[k, i] = d_k beta^i
    dbeta[2] = GEN @ beta / y
    dalpha = np.array([v['dlapse'][0], v['dlapse'][1], 0.])
    return gamma, dgamma, K, beta, dbeta, dalpha


def cartoon_triad(gamma, dgamma, gd):
    """E (frame vectors as columns), its x, y, z derivatives and time derivative, all on the plane"""
    g2 = gamma[:2, :2]
    g2i = np.linalg.inv(g2)
    dg2i = [-g2i @ dgamma[k, :2, :2] @ g2i for k in range(2)] + [-g2i @ gd[:2, :2] @ g2i]
    S, dS = sqrt2_and_deriv(g2i, dg2i)
    E = np.zeros((3, 3))
    E[:2, :2] = S
    E[2, 2] = gamma[2, 2]**-0.5
    dE = np.zeros((3, 3, 3))
    for k in range(2):
        dE[k, :2, :2] = dS[k]
        dE[k, 2, 2] = -0.5 * gamma[2, 2]**-1.5 * dgamma[k, 2, 2]
    Edot = np.zeros((3, 3))
    Edot[:2, :2] = dS[2]
    Edot[2, 2] = -0.5 * gamma[2, 2]**-1.5 * gd[2, 2]
    return E, dE, Edot


def cartoon_rhs(mu, v, y):
    gamma, dgamma, K, beta, dbeta, dalpha = cartoon_fields(v, y)
    alpha = v['lapse']
    gd = gamma_dot(alpha, beta, dbeta, gamma, dgamma, K)
    E, dE, Edot = cartoon_triad(gamma, dgamma, gd)
    dE[2] = (GEN @ E - E @ GEN) / y
    Gk, Gn = connections(alpha, dalpha, beta, dbeta, gamma, dgamma, K, E, dE, Edot)
    psi = v['psi']
    dpsi = [v['dpsi'][0], v['dpsi'][1], (I * M_J * psi + SPIN_GEN @ psi) / y]
    rhs = sum(beta[k] * dpsi[k] for k in range(2)) - alpha * Gn @ psi + alpha * mu * GAM[0] @ psi
    for b in range(3):
        nb = sum(E[k, b] * (dpsi[k] + Gk[k] @ psi) for k in range(3))
        rhs += alpha * gg(0, b + 1) @ nb
    th = np.linalg.inv(E)
    return rhs, dict(E=E, th=th, Gk=Gk, Gn=Gn, dpsi=dpsi, beta=beta)


def cartoon_sources(mu, v, y, psi_t=None):
    """rho, S_x, S_y, S_xx, S_xy, S_yy, S_ww for BOTH spinors (m = +1/2 evolved, m = -1/2 its mirror image in z)"""
    rhs, aux = cartoon_rhs(mu, v, y)
    if psi_t is None:
        psi_t = rhs
    T = stress_frame(v['psi'], frame_derivs(v['lapse'], aux['beta'], v['psi'], aux['dpsi'], psi_t, aux['E'], aux['Gk'], aux['Gn']))
    th = aux['th']
    rho = 2 * T[0, 0].real
    S_i = (-2 * th.T @ T[1:, 0].real)[:2]
    S_ij = (2 * th.T @ T[1:, 1:].real @ th)
    return rho, S_i, S_ij[:2, :2], S_ij[2, 2], 2 * np.real(v['psi'].conj() @ v['psi'])
