"""
Reference implementation (numpy) of the Proca equations on the cartoon plane (z = 0, x = symmetry axis, y = cylindrical
radius), as implemented in Source/Matter/ProcaFieldCartoon.hpp. Same equations and conventions as proca3d_ref.py.

For a twist-free axisymmetric field the azimuthal components vanish on the plane (X_z = E^z = 0, gamma_xz = gamma_yz = 0,
beta^z = 0), so the variables are X_x, X_y, E^x, E^y, X_phi and Z. The z-derivatives follow from axisymmetry
(Lie derivative along xi = y d_z - z d_y vanishes):
   d_z V^z = V^y / y for vectors (V^x, V^y and the covector components have vanishing d_z on the plane),
   d_z F^{xz} = F^{xy} / y, d_z F^{yz} = 0 for the antisymmetric F^{ij},
and the only nonzero field-strength component on the plane is F_xy (B points along z).
"""
import numpy as np


def proca_cartoon_rhs_and_sources(mu, kappa, y, alpha, dalpha, beta, dbeta, g2, dg2, gww, dgww, K2, Kww, V, dV, ddX):
    """In-plane quantities, all 2D: alpha, dalpha[k]; beta[i], dbeta[k, i] = d_k beta^i; g2[i, j] (xy block),
    dg2[k, i, j]; gww = gamma_zz, dgww[k]; K2 (xy block), Kww = K_zz. V = dict(X (2,), E (2,), Xphi, Z), complex;
    dV = dict(X[k, i], E[k, i], Xphi[k], Z[k]); ddX[k, l, i] = d_k d_l X_i. k, l, i run over x, y."""
    gi = np.linalg.inv(g2)
    dgi = np.array([-gi @ dg2[k] @ gi for k in range(2)])
    dlnsg = np.array([0.5 * np.trace(gi @ dg2[k]) + 0.5 * dgww[k] / gww for k in range(2)])
    trK = np.trace(gi @ K2) + Kww / gww
    X, E, P, Z = V['X'], V['E'], V['Xphi'], V['Z']
    dX, dE, dP, dZ = dV['X'], dV['E'], dV['Xphi'], dV['Z']

    X_up = gi @ X
    dX_up = np.array([dgi[k] @ X + gi @ dX[k] for k in range(2)])
    divX = np.trace(dX_up) + X_up @ dlnsg + X_up[1] / y
    divE = np.trace(dE) + E @ dlnsg + E[1] / y
    E_dn = g2 @ E

    q = 1. / np.linalg.det(g2)                      # F^{xy} = q F_xy
    dq = np.array([-q * np.trace(gi @ dg2[k]) for k in range(2)])
    Fxy = dX[0, 1] - dX[1, 0]
    dFxy = np.array([ddX[k, 0, 1] - ddX[k, 1, 0] for k in range(2)])
    Fup = q * Fxy
    dFup = dq * Fxy + q * dFxy
    # (1/sqrt(gamma)) d_j (sqrt(gamma) alpha F^{xy})
    dG = alpha * Fup * dlnsg + Fup * dalpha + alpha * dFup
    curl = np.array([dG[1] + alpha * Fup / y, -dG[0]])

    rhs = {}
    rhs['X'] = beta @ dX + dbeta @ X - alpha * (E_dn + dP) - P * dalpha
    rhs['E'] = beta @ dE - E @ dbeta + alpha * trK * E + alpha * mu**2 * X_up + alpha * gi @ dZ + curl
    rhs['Xphi'] = beta @ dP + alpha * trK * P - alpha * divX - X_up @ dalpha - alpha * Z
    rhs['Z'] = beta @ dZ + alpha * (divE + mu**2 * P - kappa * Z)

    re = np.real
    E2 = re(E_dn @ np.conj(E))
    B2 = q * abs(Fxy)**2
    X2 = re(X @ np.conj(X_up))
    P2 = abs(P)**2
    iso = 0.5 * E2 - 0.5 * B2 - 0.5 * mu**2 * (X2 - P2)
    FF = abs(Fxy)**2 * np.array([[gi[1, 1], -gi[0, 1]], [-gi[0, 1], gi[0, 0]]])
    src = {}
    src['rho'] = 0.5 * (E2 + B2) + 0.5 * mu**2 * (P2 + X2)
    src['Si'] = re(np.array([Fxy * np.conj(E[1]), -Fxy * np.conj(E[0])])) + mu**2 * re(P * np.conj(X))
    S = re(-np.outer(E_dn, np.conj(E_dn))) + FF + mu**2 * re(np.outer(X, np.conj(X))) + g2 * iso
    src['Sij'] = 0.5 * (S + S.T)
    src['Sww'] = gww * iso
    src['Gauss'] = divE + mu**2 * P
    src['n'] = np.imag(E @ np.conj(X))
    return rhs, src
