/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIRACFIELDCARTOON_HPP_
#define DIRACFIELDCARTOON_HPP_

#include "DimensionDefinitions.hpp"
#include "VarsTools.hpp"
#include <array>
#include <cmath>
#include <complex>

/*
 Two Dirac spinors of mass mu in an axisymmetric (twist-free) spacetime, evolved with the modified cartoon method on the
 2D grid (x = symmetry axis, y = cylindrical radius, the cartoon plane is z = 0 with z = w). This is the axisymmetric
 extension of the spherical Dirac stars of Spherical_BS_Evolver (DiracStar / ComplexSpinorField), with the same
 conventions, signature (-,+,+,+):

   (gamma^a nabla_a + mu) psi = 0,  {gamma^a, gamma^b} = 2 eta^{ab},  gamma^a = -i gamma_D^a,  psibar = psi^dag gamma_D^0,
   T_ab = (1/4)(X_ab + X_ba),  X_ab = psibar gamma_a nabla_b psi - (nabla_a psi)^dag gamma_D^0 gamma_b psi,

 with gamma_D the standard Dirac representation built from the Pauli matrices relabelled so that the spin is quantised
 along the symmetry axis: (sigma_x, sigma_y, sigma_z) = (sigma_3, sigma_1, sigma_2).

 Matter content. As for the spherical stars, the two spinors have J_x = +1/2 and J_x = -1/2. Only the first, psi, is
 evolved: the second is its mirror image in the plane z = 0 (which is a symmetry of any twist-free axisymmetric
 spacetime), so its stress tensor on the cartoon plane equals that of psi with the z-odd components (S_z, S_xz, S_yz)
 reversed. The sources are therefore 2 x the z-even part of T[psi], and the z-odd parts, which cancel, are never needed.

 Variables (4 complex = 8 real, frame components of psi on the plane):
   psi0 = upper, spin up along x    (the "g" component of the spherical star; even in y)
   psi1 = upper, spin down          (odd in y)
   psi2 = lower, spin up            (even in y)
   psi3 = lower, spin down          (odd in y)
 For the spherical star, psi = e^{-i omega t} ( g, 0, -i f (x/r), -i f (y/r) ) / sqrt(4 pi), with f, g the physical radial
 spinor functions of the 1D solver.

 Tetrad. Time gauge e_0 = n and the symmetric triad e_a^i = (gamma^{-1/2})^{ia} (Cartesian components), an algebraic
 function of the 3-metric that commutes with rotations, so it is smooth on the axis and needs no extra variables. On the
 plane gamma_ij is block diagonal (xy block and gamma_zz = h_ww/chi) and the square root of the 2x2 block has the closed
 form sqrt(M) = (M + s I)/t, s = sqrt(det M), t = sqrt(tr M + 2 s), which is differentiated analytically.

 3+1 equation (Gamma_mu = (1/8) omega_{mu ab} [gamma^a, gamma^b], K_ij = -nabla_i n_j as in GRChombo):
   d_t psi = beta^k d_k psi - alpha Gamma_n psi + alpha gamma^0 gamma^b e_b^k (d_k psi + Gamma_k psi) + alpha mu gamma^0 psi,
   Gamma_k = (1/2) omega_{k0b} gamma^0 gamma^b + (1/4) omega_{kab} gamma^a gamma^b,   omega_{k0b} = K_kj e_b^j,
   Gamma_n = (1/2) n.omega_{0b} gamma^0 gamma^b + (1/4) n.omega_{ab} gamma^a gamma^b,  n.omega_{0b} = -e_b^j d_j alpha/alpha,
   n.omega_{ab} = [ e_a.(d_t e_b - L_beta e_b) - (a <-> b) ] / (2 alpha)    (rotation of the triad as gamma_ij evolves;
                  d_t e_b follows from d_t gamma_ij = -2 alpha K_ij + L_beta gamma_ij),
   omega_{kab} = gamma_ij e_a^i (d_k e_b^j + Gamma^j_kl e_b^l).
 Cartoon: z-derivatives on the plane follow from axisymmetry, d_z T = [G, T] / y for the 3x3 tensors (all indices are
 rotated, G = dR/dphi), d_z V = G V / y for vectors, and for the spinor psi(R x) = e^{i phi/2} S(R) psi(x), so
   d_z psi = ( i/2 - (i/2) Sigma_x ) psi / y  =  (0, i psi1, 0, i psi3) / y.
 The time derivative of psi entering T_ab is the right-hand side itself (the equation of motion), so rho, S_i, S_ij need
 only data on the slice; the advection term cancels from n^a nabla_a psi and is not needed for the sources.

 All of this is checked numerically in Examples/DiracStar_2D/derivations/check_dirac_cartoon.py (against the covariant
 4D Dirac operator, against the full 3D equations on axisymmetric data, and against the 1D equations for spherical
 data), and this header is checked against that Python reference by Examples/DiracStar_2D/tests/test_dirac_cpp.
*/

//! Spinor evolution variables (enum entries c_psi0_Re ... c_psi3_Im must be defined in UserVariables.hpp)
template <class data_t> struct DiracSpinorVars
{
    data_t psi_Re[4];
    data_t psi_Im[4];

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools;
        define_enum_mapping(mapping_function, c_psi0_Re, psi_Re[0]);
        define_enum_mapping(mapping_function, c_psi0_Im, psi_Im[0]);
        define_enum_mapping(mapping_function, c_psi1_Re, psi_Re[1]);
        define_enum_mapping(mapping_function, c_psi1_Im, psi_Im[1]);
        define_enum_mapping(mapping_function, c_psi2_Re, psi_Re[2]);
        define_enum_mapping(mapping_function, c_psi2_Im, psi_Im[2]);
        define_enum_mapping(mapping_function, c_psi3_Re, psi_Re[3]);
        define_enum_mapping(mapping_function, c_psi3_Im, psi_Im[3]);
    }
};

//! The 3+1 sources of both spinors on the cartoon plane (GRChombo conventions: S_i = -gamma_i^a n^b T_ab, lower indices)
template <class data_t> struct DiracSources_t
{
    data_t rho;      //!< n^a n^b T_ab
    data_t Si[2];    //!< S_x, S_y
    data_t Sij[2][2];//!< S_xx, S_xy, S_yy
    data_t Sww;      //!< S_ww = S_zz on the plane
    data_t n;        //!< number density psi^dag psi summed over both spinors (charge = int n sqrt(gamma) d^3x)
};

namespace DiracCartoon
{
//! minimal complex arithmetic on data_t (works for doubles and GRChombo's simd types)
template <class data_t> struct cplx
{
    data_t re, im;
};
template <class data_t> inline cplx<data_t> operator+(const cplx<data_t> &a, const cplx<data_t> &b)
{
    return {a.re + b.re, a.im + b.im};
}
template <class data_t> inline cplx<data_t> operator-(const cplx<data_t> &a, const cplx<data_t> &b)
{
    return {a.re - b.re, a.im - b.im};
}
template <class data_t> inline cplx<data_t> operator*(const data_t &c, const cplx<data_t> &a)
{
    return {c * a.re, c * a.im};
}

template <class data_t> struct spinor
{
    cplx<data_t> c[4];
};
template <class data_t> inline spinor<data_t> operator+(const spinor<data_t> &a, const spinor<data_t> &b)
{
    spinor<data_t> out;
    for (int i = 0; i < 4; ++i)
        out.c[i] = a.c[i] + b.c[i];
    return out;
}
template <class data_t> inline spinor<data_t> operator*(const data_t &s, const spinor<data_t> &a)
{
    spinor<data_t> out;
    for (int i = 0; i < 4; ++i)
        out.c[i] = s * a.c[i];
    return out;
}
template <class data_t> inline spinor<data_t> zero_spinor()
{
    spinor<data_t> out;
    for (int i = 0; i < 4; ++i)
        out.c[i] = {0., 0.};
    return out;
}

using cmat = std::array<std::array<std::complex<double>, 4>, 4>;

//! (a constant complex matrix) x spinor
template <class data_t> inline spinor<data_t> apply(const cmat &M, const spinor<data_t> &p)
{
    spinor<data_t> out = zero_spinor<data_t>();
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
            const double mr = M[i][j].real(), mi = M[i][j].imag();
            if (mr != 0.)
            {
                out.c[i].re += mr * p.c[j].re;
                out.c[i].im += mr * p.c[j].im;
            }
            if (mi != 0.)
            {
                out.c[i].re -= mi * p.c[j].im;
                out.c[i].im += mi * p.c[j].re;
            }
        }
    return out;
}

//! Re < p, M q > = Re sum_ij conj(p_i) M_ij q_j
template <class data_t> inline data_t re_inner(const spinor<data_t> &p, const cmat &M, const spinor<data_t> &q)
{
    const spinor<data_t> Mq = apply(M, q);
    data_t out = 0.;
    for (int i = 0; i < 4; ++i)
        out += p.c[i].re * Mq.c[i].re + p.c[i].im * Mq.c[i].im;
    return out;
}

//! the constant matrices of the Dirac algebra in the conventions above
struct Constants
{
    cmat G0;          //!< gamma^0
    cmat G0G[3];      //!< gamma^0 gamma^b, b = x, y, z
    cmat GG[3][3];    //!< gamma^a gamma^b
    cmat BG[4];       //!< gamma_D^0 gamma_a (lower frame index), a = 0..3
    cmat spin_gen;    //!< generator of the spinor rotation about x: -(1/2) gamma^y gamma^z = -(i/2) Sigma_x

    Constants()
    {
        using cd = std::complex<double>;
        const cd I(0., 1.);
        typedef std::array<std::array<cd, 2>, 2> m2;
        const m2 s1 = {{{0., 1.}, {1., 0.}}}, s2 = {{{0., -I}, {I, 0.}}}, s3 = {{{1., 0.}, {0., -1.}}};
        const m2 sig[3] = {s3, s1, s2}; // spin quantised along the symmetry axis x
        cmat gD[4], gam[4];
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    gD[a][i][j] = 0.;
        for (int i = 0; i < 2; ++i)
        {
            gD[0][i][i] = 1.;
            gD[0][i + 2][i + 2] = -1.;
        }
        for (int k = 0; k < 3; ++k)
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                {
                    gD[k + 1][i][j + 2] = sig[k][i][j];
                    gD[k + 1][i + 2][j] = -sig[k][i][j];
                }
        for (int a = 0; a < 4; ++a)
            gam[a] = scale(-I, gD[a]);
        G0 = gam[0];
        for (int b = 0; b < 3; ++b)
        {
            G0G[b] = mult(gam[0], gam[b + 1]);
            for (int a = 0; a < 3; ++a)
                GG[a][b] = mult(gam[a + 1], gam[b + 1]);
        }
        BG[0] = mult(gD[0], scale(-1., gam[0]));
        for (int a = 1; a < 4; ++a)
            BG[a] = mult(gD[0], gam[a]);
        spin_gen = scale(-0.5, mult(gam[2], gam[3]));
    }

    static cmat mult(const cmat &A, const cmat &B)
    {
        cmat C;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            {
                C[i][j] = 0.;
                for (int k = 0; k < 4; ++k)
                    C[i][j] += A[i][k] * B[k][j];
            }
        return C;
    }
    static cmat scale(std::complex<double> c, const cmat &A)
    {
        cmat C;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                C[i][j] = c * A[i][j];
        return C;
    }
};

inline const Constants &constants()
{
    static const Constants c;
    return c;
}
} // namespace DiracCartoon

class DiracFieldCartoon
{
  public:
    struct params_t
    {
        double mass; //!< fermion mass mu
    };

    explicit DiracFieldCartoon(params_t a_params) : m_params(a_params) {}

    //! Right-hand side of the spinor (excluding the advection term beta^k d_k psi, which the caller adds from its upwinded
    //! advection derivatives) and the sources of both spinors. vars_t must hold the CCZ4 cartoon variables (chi, h, hww,
    //! K, A, Aww, lapse, shift) and psi_Re[4], psi_Im[4]; d1 their first derivatives; y the cartoon coordinate (> 0).
    template <class data_t, template <typename> class vars_t>
    void compute(DiracSpinorVars<data_t> &rhs_no_advec, DiracSources_t<data_t> &src, const vars_t<data_t> &vars,
                 const vars_t<Tensor<1, data_t>> &d1, const double y) const;

    //! number density of both spinors, 2 psi^dag psi
    template <class data_t, template <typename> class vars_t>
    static data_t number_density(const vars_t<data_t> &vars)
    {
        data_t n = 0.;
        for (int c = 0; c < 4; ++c)
            n += vars.psi_Re[c] * vars.psi_Re[c] + vars.psi_Im[c] * vars.psi_Im[c];
        return 2. * n;
    }

  protected:
    params_t m_params;
};

template <class data_t, template <typename> class vars_t>
void DiracFieldCartoon::compute(DiracSpinorVars<data_t> &rhs_no_advec, DiracSources_t<data_t> &src,
                                const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
                                const double y) const
{
    using namespace DiracCartoon;
    const Constants &C = constants();
    const double GEN[3][3] = {{0., 0., 0.}, {0., 0., -1.}, {0., 1., 0.}}; // d R_phi / d phi, rotation about x
    const double oy = 1. / y;
    const data_t alpha = vars.lapse;
    const data_t chi = vars.chi;
    const double mu = m_params.mass;

    // ---------------------------------------------------------------- 3x3 metric data on the plane
    data_t g[3][3], dg[3][3][3], K[3][3], gd[3][3], gi[3][3];
    data_t beta[3], dbeta[3][3], dalpha[3]; // dbeta[k][i] = d_k beta^i
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            g[i][j] = 0.;
            K[i][j] = 0.;
            for (int k = 0; k < 3; ++k)
                dg[k][i][j] = 0.;
        }
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
            g[i][j] = vars.h[i][j] / chi;
            K[i][j] = (vars.A[i][j] + vars.h[i][j] * vars.K / 3.) / chi;
            for (int k = 0; k < 2; ++k)
                dg[k][i][j] = (d1.h[i][j][k] - vars.h[i][j] * d1.chi[k] / chi) / chi;
        }
    g[2][2] = vars.hww / chi;
    K[2][2] = (vars.Aww + vars.hww * vars.K / 3.) / chi;
    for (int k = 0; k < 2; ++k)
        dg[k][2][2] = (d1.hww[k] - vars.hww * d1.chi[k] / chi) / chi;
    // d_z gamma = [G, gamma] / y
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            data_t c = 0.;
            for (int l = 0; l < 3; ++l)
                c += GEN[i][l] * g[l][j] - g[i][l] * GEN[l][j];
            dg[2][i][j] = c * oy;
        }

    beta[0] = vars.shift[0];
    beta[1] = vars.shift[1];
    beta[2] = 0.;
    for (int k = 0; k < 3; ++k)
        for (int i = 0; i < 3; ++i)
            dbeta[k][i] = 0.;
    for (int k = 0; k < 2; ++k)
        for (int i = 0; i < 2; ++i)
            dbeta[k][i] = d1.shift[i][k];
    for (int i = 0; i < 3; ++i)
    {
        data_t c = 0.;
        for (int l = 0; l < 3; ++l)
            c += GEN[i][l] * beta[l];
        dbeta[2][i] = c * oy;
    }
    dalpha[0] = d1.lapse[0];
    dalpha[1] = d1.lapse[1];
    dalpha[2] = 0.;

    // d_t gamma_ij = -2 alpha K_ij + beta^k d_k gamma_ij + gamma_kj d_i beta^k + gamma_ik d_j beta^k
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            gd[i][j] = -2. * alpha * K[i][j];
            for (int k = 0; k < 3; ++k)
                gd[i][j] += beta[k] * dg[k][i][j] + g[k][j] * dbeta[i][k] + g[i][k] * dbeta[j][k];
        }

    // inverse metric (block diagonal)
    const data_t det2 = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            gi[i][j] = 0.;
    gi[0][0] = g[1][1] / det2;
    gi[1][1] = g[0][0] / det2;
    gi[0][1] = -g[0][1] / det2;
    gi[1][0] = gi[0][1];
    gi[2][2] = 1. / g[2][2];

    // Christoffel symbols of gamma, G3[i][j][k] = Gamma^i_jk
    data_t G3[3][3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
            {
                data_t c = 0.;
                for (int l = 0; l < 3; ++l)
                    c += gi[i][l] * (dg[j][l][k] + dg[k][l][j] - dg[l][j][k]);
                G3[i][j][k] = 0.5 * c;
            }

    // ---------------------------------------------------------------- symmetric triad and its derivatives
    // M = (gamma^{-1}) xy block, E2 = sqrt(M) = (M + s I)/t; derivative directions x, y and t
    data_t E[3][3], dE[3][3][3], Edot[3][3], th[3][3];
    {
        const data_t s = sqrt(gi[0][0] * gi[1][1] - gi[0][1] * gi[1][0]);
        const data_t t = sqrt(gi[0][0] + gi[1][1] + 2. * s);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                E[i][j] = 0.;
                Edot[i][j] = 0.;
                for (int k = 0; k < 3; ++k)
                    dE[k][i][j] = 0.;
            }
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                E[i][j] = (gi[i][j] + ((i == j) ? s : data_t(0.))) / t;

        // dM = -M dgamma2 M for each direction (dir 2 = time)
        for (int dir = 0; dir < 3; ++dir)
        {
            data_t dg2[2][2], dM[2][2];
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    dg2[i][j] = (dir < 2) ? dg[dir][i][j] : gd[i][j];
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                {
                    data_t c = 0.;
                    for (int k = 0; k < 2; ++k)
                        for (int l = 0; l < 2; ++l)
                            c -= gi[i][k] * dg2[k][l] * gi[l][j];
                    dM[i][j] = c;
                }
            // ds = (s/2) tr(M^{-1} dM), M^{-1} = gamma2
            data_t trMinvdM = 0.;
            for (int i = 0; i < 2; ++i)
                for (int k = 0; k < 2; ++k)
                    trMinvdM += g[i][k] * dM[k][i];
            const data_t ds = 0.5 * s * trMinvdM;
            const data_t dt = (dM[0][0] + dM[1][1] + 2. * ds) / (2. * t);
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                {
                    const data_t val = (dM[i][j] + ((i == j) ? ds : data_t(0.))) / t - E[i][j] * dt / t;
                    if (dir < 2)
                        dE[dir][i][j] = val;
                    else
                        Edot[i][j] = val;
                }
        }
        E[2][2] = 1. / sqrt(g[2][2]);
        const data_t m32 = -0.5 * E[2][2] / g[2][2]; // d(g^{-1/2})/dg = -1/2 g^{-3/2}
        dE[0][2][2] = m32 * dg[0][2][2];
        dE[1][2][2] = m32 * dg[1][2][2];
        Edot[2][2] = m32 * gd[2][2];
        // d_z E = [G, E] / y
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                data_t c = 0.;
                for (int l = 0; l < 3; ++l)
                    c += GEN[i][l] * E[l][j] - E[i][l] * GEN[l][j];
                dE[2][i][j] = c * oy;
            }
        // co-triad theta = E^{-1} = gamma^{1/2}
        const data_t detE = E[0][0] * E[1][1] - E[0][1] * E[1][0];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                th[i][j] = 0.;
        th[0][0] = E[1][1] / detE;
        th[1][1] = E[0][0] / detE;
        th[0][1] = -E[0][1] / detE;
        th[1][0] = -E[1][0] / detE;
        th[2][2] = 1. / E[2][2];
    }

    // ---------------------------------------------------------------- spin connection
    data_t om[3][3][3], om0[3][3], an0[3], anab[3][3];
    for (int k = 0; k < 3; ++k)
        for (int b = 0; b < 3; ++b)
        {
            // (nabla_k e_b)^j = d_k e_b^j + Gamma^j_kl e_b^l
            data_t nab[3];
            for (int j = 0; j < 3; ++j)
            {
                nab[j] = dE[k][j][b];
                for (int l = 0; l < 3; ++l)
                    nab[j] += G3[j][k][l] * E[l][b];
            }
            for (int a = 0; a < 3; ++a)
            {
                data_t c = 0.;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        c += g[i][j] * E[i][a] * nab[j];
                om[k][a][b] = c;
            }
            data_t c0 = 0.;
            for (int j = 0; j < 3; ++j)
                c0 += K[k][j] * E[j][b];
            om0[k][b] = c0;
        }
    for (int b = 0; b < 3; ++b)
    {
        data_t c = 0.;
        for (int j = 0; j < 3; ++j)
            c -= E[j][b] * dalpha[j];
        an0[b] = c / alpha;
    }
    {
        // W_ab = e_a . (d_t e_b - L_beta e_b),  (L_beta e_b)^i = beta^l d_l e_b^i - e_b^l d_l beta^i
        data_t V[3][3]; // V[i][b] = d_t e_b^i - (L_beta e_b)^i
        for (int i = 0; i < 3; ++i)
            for (int b = 0; b < 3; ++b)
            {
                data_t c = Edot[i][b];
                for (int l = 0; l < 3; ++l)
                    c -= beta[l] * dE[l][i][b] - E[l][b] * dbeta[l][i];
                V[i][b] = c;
            }
        data_t W[3][3];
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
            {
                data_t c = 0.;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        c += E[i][a] * g[i][j] * V[j][b];
                W[a][b] = c;
            }
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                anab[a][b] = (W[a][b] - W[b][a]) / (2. * alpha);
    }

    // ---------------------------------------------------------------- spinor and its derivatives
    spinor<data_t> psi, dpsi[3];
    for (int c = 0; c < 4; ++c)
    {
        psi.c[c] = {vars.psi_Re[c], vars.psi_Im[c]};
        for (int k = 0; k < 2; ++k)
            dpsi[k].c[c] = {d1.psi_Re[c][k], d1.psi_Im[c][k]};
    }
    // d_z psi = (i m + spin_gen) psi / y with m = +1/2
    {
        spinor<data_t> sg = apply(C.spin_gen, psi);
        for (int c = 0; c < 4; ++c)
            dpsi[2].c[c] = {(sg.c[c].re - 0.5 * psi.c[c].im) * oy, (sg.c[c].im + 0.5 * psi.c[c].re) * oy};
    }

    // Gamma_k psi and Gamma_n psi
    spinor<data_t> Gk_psi[3], Gn_psi = zero_spinor<data_t>();
    {
        spinor<data_t> G0Gpsi[3], GGpsi[3][3];
        for (int b = 0; b < 3; ++b)
        {
            G0Gpsi[b] = apply(C.G0G[b], psi);
            for (int a = 0; a < 3; ++a)
                GGpsi[a][b] = (a == b) ? psi : apply(C.GG[a][b], psi);
        }
        for (int k = 0; k < 3; ++k)
        {
            Gk_psi[k] = zero_spinor<data_t>();
            for (int b = 0; b < 3; ++b)
            {
                Gk_psi[k] = Gk_psi[k] + data_t(0.5 * om0[k][b]) * G0Gpsi[b];
                for (int a = 0; a < 3; ++a)
                    if (a != b)
                        Gk_psi[k] = Gk_psi[k] + data_t(0.25 * om[k][a][b]) * GGpsi[a][b];
            }
        }
        for (int b = 0; b < 3; ++b)
        {
            Gn_psi = Gn_psi + data_t(0.5 * an0[b]) * G0Gpsi[b];
            for (int a = 0; a < 3; ++a)
                if (a != b)
                    Gn_psi = Gn_psi + data_t(0.25 * anab[a][b]) * GGpsi[a][b];
        }
    }

    // frame derivatives along the spatial triad, nb[b] = e_b^k (d_k psi + Gamma_k psi)
    spinor<data_t> nb[3];
    for (int b = 0; b < 3; ++b)
    {
        nb[b] = zero_spinor<data_t>();
        for (int k = 0; k < 3; ++k)
            nb[b] = nb[b] + E[k][b] * (dpsi[k] + Gk_psi[k]);
    }

    // ---------------------------------------------------------------- right-hand side (without advection)
    spinor<data_t> rhs = (-1. * alpha) * Gn_psi + (alpha * mu) * apply(C.G0, psi);
    for (int b = 0; b < 3; ++b)
        rhs = rhs + alpha * apply(C.G0G[b], nb[b]);
    for (int c = 0; c < 4; ++c)
    {
        rhs_no_advec.psi_Re[c] = rhs.c[c].re;
        rhs_no_advec.psi_Im[c] = rhs.c[c].im;
    }

    // ---------------------------------------------------------------- sources
    // nabla_0 psi = n^mu nabla_mu psi = (d_t psi - beta^k d_k psi)/alpha + Gamma_n psi
    spinor<data_t> nab[4];
    nab[0] = (1. / alpha) * rhs + Gn_psi;
    for (int b = 0; b < 3; ++b)
        nab[b + 1] = nb[b];
    // T_ab = (1/4) Re (X_ab + X_ba),  Re X_ab = Re <psi, BG_a nab_b> - Re <nab_a, BG_b psi>
    auto T = [&](int a, int b) {
        return 0.25 * (re_inner(psi, C.BG[a], nab[b]) - re_inner(nab[a], C.BG[b], psi) +
                       re_inner(psi, C.BG[b], nab[a]) - re_inner(nab[b], C.BG[a], psi));
    };
    // both spinors: 2 x (z-even part of T[psi])
    src.rho = 2. * T(0, 0);
    data_t T0[2] = {T(1, 0), T(2, 0)};
    data_t Tab[2][2] = {{T(1, 1), T(1, 2)}, {T(1, 2), T(2, 2)}};
    for (int i = 0; i < 2; ++i)
    {
        src.Si[i] = -2. * (th[0][i] * T0[0] + th[1][i] * T0[1]);
        for (int j = 0; j < 2; ++j)
        {
            data_t c = 0.;
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    c += th[a][i] * th[b][j] * Tab[a][b];
            src.Sij[i][j] = 2. * c;
        }
    }
    src.Sww = 2. * th[2][2] * th[2][2] * T(3, 3);
    src.n = number_density(vars);
}

#endif /* DIRACFIELDCARTOON_HPP_ */
