/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCAFIELDCARTOON_HPP_
#define PROCAFIELDCARTOON_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include <cmath>

/*
 A complex Proca field of mass mu in an axisymmetric (twist-free) spacetime, evolved with the modified cartoon method on
 the 2D grid (x = symmetry axis, y = cylindrical radius, the cartoon plane is z = 0 with z = w). This is the
 axisymmetric extension of the spherical Proca stars of Spherical_BS_Evolver (ProcaStar / ComplexVectorField, branch
 proca-stars), with the same conventions, signature (-,+,+,+), K_ij = -nabla_i n_j:

   L = -1/4 F_ab Fbar^ab - 1/2 mu^2 A_a Abar^a,   nabla_a F^{ab} = mu^2 A^b.

 Variables (Zilhao, Witek & Cardoso 2015), real and imaginary parts of each: X_phi = -n^a A_a, X_i = gamma_i^a A_a,
 E^i = gamma^{ij} F_{ja} n^a and the Gauss-constraint damping variable Z. For a twist-free axisymmetric field the
 azimuthal components vanish on the plane (X_z = E^z = 0), so the evolved components are X_x, X_y, E^x, E^y, X_phi, Z
 (12 real variables). Equations (F_ij = d_i X_j - d_j X_i):
   d_t X_i   = L_beta X_i - alpha (E_i + D_i X_phi) - X_phi D_i alpha
   d_t E^i   = L_beta E^i + alpha K E^i + D_j(alpha F^{ij}) + alpha mu^2 X^i + alpha D^i Z
   d_t X_phi = L_beta X_phi + alpha K X_phi - alpha D_i X^i - X^i D_i alpha - alpha Z
   d_t Z     = L_beta Z + alpha (D_i E^i + mu^2 X_phi - kappa Z)
 Sources (sums over the real and imaginary parts), B^2 = F_ij F^ij / 2:
   rho  = (E_i E^i + B^2)/2 + mu^2 (X_phi^2 + X_i X^i)/2,     S_i = F_ij E^j + mu^2 X_phi X_i,
   S_ij = -E_i E_j + F_ik F_jl gamma^kl + mu^2 X_i X_j + gamma_ij [E^2/2 - B^2/2 - mu^2 (X_k X^k - X_phi^2)/2].
 Number density n = Im(E^i Xbar_i) (normal observers; particle number int n sqrt(gamma) d^3x), Gauss constraint
 C = D_i E^i + mu^2 X_phi.

 Cartoon terms (from the vanishing Lie derivative along xi = y d_z - z d_y): on the plane the only nonzero component
 of F_ij is F_xy (the "magnetic field" points along z), F^{xy} = F_xy / det(gamma_xy block),
   D_i V^i = (1/sqrt(gamma)) d_a (sqrt(gamma) V^a) + V^y / y            (a = x, y; d_z V^z = V^y / y),
   D_j(alpha F^{xj}) = (1/sqrt(gamma)) d_y(sqrt(gamma) alpha F^{xy}) + alpha F^{xy} / y   (d_z F^{xz} = F^{xy} / y),
   D_j(alpha F^{yj}) = -(1/sqrt(gamma)) d_x(sqrt(gamma) alpha F^{xy}),
 sqrt(gamma) = sqrt(det(gamma_xy block) gamma_ww). The Lie derivatives of the in-plane components need no cartoon
 terms, and S_ww = gamma_ww [E^2/2 - B^2/2 - mu^2 (X_k X^k - X_phi^2)/2].

 Only the free field is evolved, as in the 1D code. Checked numerically in Examples/ProcaStar_2D/derivations
 (check_proca_cartoon.py against the full 3D equations on axisymmetric data; check_proca3d.py for the 3D equations
 against the covariant 4D theory and the 1D code), and this header against the Python reference by
 Examples/ProcaStar_2D/tests/test_proca_cpp.
*/

//! Proca evolution variables (enum entries c_Xx_Re ... c_Z_Im must be defined in UserVariables.hpp)
template <class data_t> struct ProcaCartoonVars
{
    data_t X_Re[2], X_Im[2]; //!< X_x, X_y
    data_t E_Re[2], E_Im[2]; //!< E^x, E^y
    data_t Xphi_Re, Xphi_Im, Z_Re, Z_Im;

    template <typename mapping_function_t> void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools;
        define_enum_mapping(mapping_function, c_Xx_Re, X_Re[0]);
        define_enum_mapping(mapping_function, c_Xy_Re, X_Re[1]);
        define_enum_mapping(mapping_function, c_Xx_Im, X_Im[0]);
        define_enum_mapping(mapping_function, c_Xy_Im, X_Im[1]);
        define_enum_mapping(mapping_function, c_Ex_Re, E_Re[0]);
        define_enum_mapping(mapping_function, c_Ey_Re, E_Re[1]);
        define_enum_mapping(mapping_function, c_Ex_Im, E_Im[0]);
        define_enum_mapping(mapping_function, c_Ey_Im, E_Im[1]);
        define_enum_mapping(mapping_function, c_Xphi_Re, Xphi_Re);
        define_enum_mapping(mapping_function, c_Xphi_Im, Xphi_Im);
        define_enum_mapping(mapping_function, c_Z_Re, Z_Re);
        define_enum_mapping(mapping_function, c_Z_Im, Z_Im);
    }
};

//! the Proca variables needing second derivatives (X_i, for the curl term)
template <class data_t> struct ProcaCartoonDiff2Vars
{
    data_t X_Re[2], X_Im[2];

    template <typename mapping_function_t> void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools;
        define_enum_mapping(mapping_function, c_Xx_Re, X_Re[0]);
        define_enum_mapping(mapping_function, c_Xy_Re, X_Re[1]);
        define_enum_mapping(mapping_function, c_Xx_Im, X_Im[0]);
        define_enum_mapping(mapping_function, c_Xy_Im, X_Im[1]);
    }
};

//! The 3+1 sources on the cartoon plane (GRChombo conventions: S_i = -gamma_i^a n^b T_ab, lower indices) and diagnostics
template <class data_t> struct ProcaSources_t
{
    data_t rho;       //!< n^a n^b T_ab
    data_t Si[2];     //!< S_x, S_y
    data_t Sij[2][2]; //!< S_xx, S_xy, S_yy
    data_t Sww;       //!< S_ww = S_zz on the plane
    data_t n;         //!< number density Im(E^i Xbar_i)
    data_t C_Re, C_Im;//!< Gauss constraint D_i E^i + mu^2 X_phi
};

class ProcaFieldCartoon
{
  public:
    struct params_t
    {
        double mass;    //!< mu
        double damping; //!< kappa
    };

    explicit ProcaFieldCartoon(params_t a_params) : m_params(a_params) {}

    //! metric data on the plane
    template <class data_t> struct metric_t
    {
        data_t g[2][2], gi[2][2], dgi[2][2][2]; //!< xy block, its inverse, d_k of the inverse ([k][i][j])
        data_t gww;                             //!< gamma_ww
        data_t dlnsg[2];                        //!< d_k ln sqrt(gamma)
        data_t q, dq[2];                        //!< 1/det(xy block) and its derivatives
    };

    template <class data_t, template <typename> class vars_t>
    static metric_t<data_t> metric(const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1)
    {
        metric_t<data_t> m;
        const data_t chi = vars.chi;
        data_t dg[2][2][2];
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
            {
                m.g[i][j] = vars.h[i][j] / chi;
                for (int k = 0; k < 2; ++k)
                    dg[k][i][j] = (d1.h[i][j][k] - vars.h[i][j] * d1.chi[k] / chi) / chi;
            }
        m.gww = vars.hww / chi;
        const data_t det2 = m.g[0][0] * m.g[1][1] - m.g[0][1] * m.g[1][0];
        m.q = 1. / det2;
        m.gi[0][0] = m.g[1][1] * m.q;
        m.gi[1][1] = m.g[0][0] * m.q;
        m.gi[0][1] = -m.g[0][1] * m.q;
        m.gi[1][0] = m.gi[0][1];
        for (int k = 0; k < 2; ++k)
        {
            const data_t dgww = (d1.hww[k] - vars.hww * d1.chi[k] / chi) / chi;
            data_t tr = 0.;
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                {
                    tr += m.gi[i][j] * dg[k][j][i];
                    data_t c = 0.;
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b)
                            c -= m.gi[i][a] * dg[k][a][b] * m.gi[b][j];
                    m.dgi[k][i][j] = c;
                }
            m.dlnsg[k] = 0.5 * tr + 0.5 * dgww / m.gww;
            m.dq[k] = -m.q * tr;
        }
        return m;
    }

    //! Sources and diagnostics (no lapse or shift needed)
    template <class data_t, template <typename> class vars_t>
    void compute_sources(ProcaSources_t<data_t> &src, const vars_t<data_t> &vars,
                         const vars_t<Tensor<1, data_t>> &d1, const double y) const
    {
        const auto m = metric(vars, d1);
        src.rho = 0.;
        src.Sww = 0.;
        for (int i = 0; i < 2; ++i)
        {
            src.Si[i] = 0.;
            for (int j = 0; j < 2; ++j)
                src.Sij[i][j] = 0.;
        }
        add_component_sources(src, m, vars.X_Re, vars.E_Re, vars.Xphi_Re, d1.X_Re);
        add_component_sources(src, m, vars.X_Im, vars.E_Im, vars.Xphi_Im, d1.X_Im);
        src.n = 0.;
        for (int i = 0; i < 2; ++i)
            src.n += vars.E_Im[i] * vars.X_Re[i] - vars.E_Re[i] * vars.X_Im[i];
        const double mu2 = m_params.mass * m_params.mass;
        src.C_Re = mu2 * vars.Xphi_Re + vars.E_Re[1] / y;
        src.C_Im = mu2 * vars.Xphi_Im + vars.E_Im[1] / y;
        for (int i = 0; i < 2; ++i)
        {
            src.C_Re += d1.E_Re[i][i] + vars.E_Re[i] * m.dlnsg[i];
            src.C_Im += d1.E_Im[i][i] + vars.E_Im[i] * m.dlnsg[i];
        }
    }

    //! Right-hand side of the Proca variables excluding the advection terms beta^k d_k (which the caller adds from its
    //! upwinded derivatives), and the sources. vars_t must hold the CCZ4 cartoon variables (chi, h, hww, K, lapse,
    //! shift) and the Proca variables, d1 their first derivatives, d2 the second derivatives of X_i; y > 0.
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    void compute(ProcaCartoonVars<data_t> &rhs, ProcaSources_t<data_t> &src, const vars_t<data_t> &vars,
                 const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2, const double y) const
    {
        compute_sources(src, vars, d1, y);
        const auto m = metric(vars, d1);
        component_rhs(rhs.X_Re, rhs.E_Re, rhs.Xphi_Re, rhs.Z_Re, vars, d1, m, vars.X_Re, vars.E_Re, vars.Xphi_Re,
                      vars.Z_Re, d1.X_Re, d1.E_Re, d1.Xphi_Re, d1.Z_Re, d2.X_Re, y);
        component_rhs(rhs.X_Im, rhs.E_Im, rhs.Xphi_Im, rhs.Z_Im, vars, d1, m, vars.X_Im, vars.E_Im, vars.Xphi_Im,
                      vars.Z_Im, d1.X_Im, d1.E_Im, d1.Xphi_Im, d1.Z_Im, d2.X_Im, y);
    }

  protected:
    params_t m_params;

    template <class data_t, class d1_t>
    void add_component_sources(ProcaSources_t<data_t> &src, const metric_t<data_t> &m, const data_t *X,
                               const data_t *E, const data_t &P, const d1_t &dX) const
    {
        const double mu2 = m_params.mass * m_params.mass;
        const data_t Fxy = dX[1][0] - dX[0][1]; // d_x X_y - d_y X_x
        data_t E_dn[2], X_up[2];
        for (int i = 0; i < 2; ++i)
        {
            E_dn[i] = m.g[i][0] * E[0] + m.g[i][1] * E[1];
            X_up[i] = m.gi[i][0] * X[0] + m.gi[i][1] * X[1];
        }
        const data_t E2 = E_dn[0] * E[0] + E_dn[1] * E[1];
        const data_t X2 = X[0] * X_up[0] + X[1] * X_up[1];
        const data_t B2 = m.q * Fxy * Fxy;
        const data_t iso = 0.5 * E2 - 0.5 * B2 - 0.5 * mu2 * (X2 - P * P);
        const data_t FF[2][2] = {{Fxy * Fxy * m.gi[1][1], -Fxy * Fxy * m.gi[0][1]},
                                 {-Fxy * Fxy * m.gi[0][1], Fxy * Fxy * m.gi[0][0]}};
        src.rho += 0.5 * (E2 + B2) + 0.5 * mu2 * (P * P + X2);
        src.Si[0] += Fxy * E[1] + mu2 * P * X[0];
        src.Si[1] += -Fxy * E[0] + mu2 * P * X[1];
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                src.Sij[i][j] += -E_dn[i] * E_dn[j] + FF[i][j] + mu2 * X[i] * X[j] + m.g[i][j] * iso;
        src.Sww += m.gww * iso;
    }

    template <class data_t, template <typename> class vars_t, class d1_t, class d2_t>
    void component_rhs(data_t *rX, data_t *rE, data_t &rP, data_t &rZ, const vars_t<data_t> &vars,
                       const vars_t<Tensor<1, data_t>> &d1, const metric_t<data_t> &m, const data_t *X, const data_t *E,
                       const data_t &P, const data_t &Z, const d1_t &dX, const d1_t &dE, const Tensor<1, data_t> &dP,
                       const Tensor<1, data_t> &dZ, const d2_t &ddX, const double y) const
    {
        const double mu2 = m_params.mass * m_params.mass;
        const double kappa = m_params.damping;
        const data_t alpha = vars.lapse;
        const data_t trK = vars.K;

        data_t E_dn[2], X_up[2];
        for (int i = 0; i < 2; ++i)
        {
            E_dn[i] = m.g[i][0] * E[0] + m.g[i][1] * E[1];
            X_up[i] = m.gi[i][0] * X[0] + m.gi[i][1] * X[1];
        }
        data_t divX = X_up[1] / y, divE = E[1] / y;
        for (int i = 0; i < 2; ++i)
        {
            divX += X_up[i] * m.dlnsg[i];
            divE += dE[i][i] + E[i] * m.dlnsg[i];
            for (int j = 0; j < 2; ++j)
                divX += m.dgi[i][i][j] * X[j] + m.gi[i][j] * dX[j][i]; // d_i X^i
        }

        // curl: G = sqrt(gamma) alpha F^{xy}, (1/sqrt(gamma)) d_k G = alpha F^{xy} dlnsg_k + F^{xy} d_k alpha + alpha d_k F^{xy}
        const data_t Fxy = dX[1][0] - dX[0][1];
        const data_t Fup = m.q * Fxy;
        data_t dG[2];
        for (int k = 0; k < 2; ++k)
        {
            const data_t dFxy = ddX[1][0][k] - ddX[0][1][k];
            dG[k] = alpha * Fup * m.dlnsg[k] + Fup * d1.lapse[k] + alpha * (m.dq[k] * Fxy + m.q * dFxy);
        }
        const data_t curl[2] = {dG[1] + alpha * Fup / y, -dG[0]};

        for (int i = 0; i < 2; ++i)
        {
            rX[i] = -alpha * (E_dn[i] + dP[i]) - P * d1.lapse[i];
            rE[i] = alpha * trK * E[i] + alpha * mu2 * X_up[i] + curl[i];
            for (int j = 0; j < 2; ++j)
            {
                rX[i] += X[j] * d1.shift[j][i];  // X_j d_i beta^j
                rE[i] += -E[j] * d1.shift[i][j]; // -E^j d_j beta^i
                rE[i] += alpha * m.gi[i][j] * dZ[j];
            }
        }
        rP = alpha * trK * P - alpha * divX - alpha * Z - X_up[0] * d1.lapse[0] - X_up[1] * d1.lapse[1];
        rZ = alpha * (divE + mu2 * P - kappa * Z);
    }
};

#endif /* PROCAFIELDCARTOON_HPP_ */
