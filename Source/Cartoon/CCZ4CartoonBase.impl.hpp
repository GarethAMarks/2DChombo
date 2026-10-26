/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4CARTOONBASE_HPP_)
#error "This file should only be included through CCZ4CartoonBase.hpp"
#endif

#ifndef CCZ4CARTOONBASE_IMPL_HPP_
#define CCZ4CARTOONBASE_IMPL_HPP_

#define COVARIANTZ4
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

template <class gauge_t, class deriv_t>
inline CCZ4CartoonBase<gauge_t, deriv_t>::CCZ4CartoonBase(params_t params, double dx, double sigma,
                                                        double a_G_Newton, int formulation,
                                                        double cosmological_constant)
    : CCZ4RHS<gauge_t, deriv_t>(params, dx, sigma, formulation, cosmological_constant), m_G_Newton(a_G_Newton)
{
}

// Vacuum CCZ4 right-hand side with the cartoon terms, including the gauge (moved verbatim from CCZ4Cartoon::rhs_equation)
template <class gauge_t, class deriv_t>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
void CCZ4CartoonBase<gauge_t, deriv_t>::vacuum_rhs_equation(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                                                           const vars_t<Tensor<1, data_t>> &d1,
                                                           const diff2_vars_t<Tensor<2, data_t>> &d2,
                                                           const vars_t<data_t> &advec,
                                                           const double &cartoon_coord) const
{
using namespace TensorAlgebra;
const int dI = CH_SPACEDIM - 1;
const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
const double two_over_gr_spacedim = 2. * one_over_gr_spacedim;

auto h_UU = compute_inverse_sym(vars.h);
auto h_UU_ww = 1. / vars.hww;
auto chris = compute_christoffel(d1.h, h_UU);
const double one_over_cartoon_coord = 1. / cartoon_coord;
const double one_over_cartoon_coord2 = 1. / (cartoon_coord * cartoon_coord);

Tensor<1, data_t> chris_ww;
FOR(i)
        {
                chris_ww[i] =
                        one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
        }
Tensor<1, data_t> chris_contracted; //!< includes the higher D contributions
FOR(i)
        {
                chris_contracted[i] = chris.contracted[i] + nS * h_UU_ww * chris_ww[i];
        }

Tensor<1, data_t> Z_over_chi;
Tensor<1, data_t> Z;
if (this->m_formulation == CCZ4::USE_BSSN)
{
FOR(i) Z_over_chi[i] = 0.0;
}
else
{
FOR(i) Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris_contracted[i]);
}
FOR(i) Z[i] = vars.chi * Z_over_chi[i];

auto ricci = CCZ4CartoonGeometry::compute_ricci_Z(
        vars, d1, d2, h_UU, h_UU_ww, chris, Z_over_chi, cartoon_coord);

data_t divshift = compute_trace(d1.shift) +
                  nS * vars.shift[dI] /
                  cartoon_coord; //!< includes the higher D contribution
data_t divshift_w = compute_trace(d1.shift) -
                    CH_SPACEDIM * vars.shift[dI] /
                    cartoon_coord; //!< combination that appears in the
//!< rhs of the ww components
data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);
data_t dlapse_dot_dhww = compute_dot_product(d1.lapse, d1.hww, h_UU);

data_t hUU_dlapse_cartoon = 0;
FOR(i)
        {
                hUU_dlapse_cartoon +=
                        one_over_cartoon_coord * h_UU[dI][i] * d1.lapse[i];
        }

Tensor<1, data_t> reg_03, reg_04;
FOR(i)
        {
                reg_03[i] = one_over_cartoon_coord * d1.shift[i][dI] -
                            one_over_cartoon_coord2 * delta(i, dI) * vars.shift[dI];
        reg_04[i] = one_over_cartoon_coord * d1.shift[dI][i] -
        one_over_cartoon_coord2 * delta(i, dI) * vars.shift[dI];
        }

Tensor<2, data_t> covdtilde2lapse;
Tensor<2, data_t> covd2lapse; // NOTE: we compute chi * D_i D_j lapse
FOR(k, l)
{
covdtilde2lapse[k][l] = d2.lapse[k][l];
FOR(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
covd2lapse[k][l] =
vars.chi * covdtilde2lapse[k][l] +
0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
vars.h[k][l] * dlapse_dot_dchi);
}
data_t covd2lapse_ww =
        vars.chi * (0.5 * dlapse_dot_dhww + vars.hww * hUU_dlapse_cartoon) -
        0.5 * vars.hww * dlapse_dot_dchi;

/*data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
FOR(i)
{
    tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
    FOR(j)
    {
        tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                       d1.lapse[i] * d1.chi[j]);
    }
}*/
data_t tr_covd2lapse = compute_trace(covd2lapse, h_UU);
tr_covd2lapse += nS * h_UU_ww * covd2lapse_ww;

// Make A_{ij} trace free
/*data_t trA = nS * h_UU_ww * vars.Aww;
FOR(i,j)
{
    trA += h_UU[i][j] * vars.A[i][j];
}*/
data_t trA = compute_trace(vars.A, h_UU);
trA += nS * h_UU_ww * vars.Aww;

Tensor<2, data_t> A_TF;
FOR(i, j)
{
A_TF[i][j] = vars.A[i][j] - trA * vars.h[i][j] * one_over_gr_spacedim;
}
data_t Aww_TF;
Aww_TF = vars.Aww - trA * vars.hww * one_over_gr_spacedim;

Tensor<2, data_t> A_UU = raise_all(A_TF, h_UU);
data_t A_UU_ww = h_UU_ww * h_UU_ww * Aww_TF;
// A^{ij} A_{ij}. - Note the abuse of the compute trace function.
data_t tr_A2 = compute_trace(A_TF, A_UU) + nS * A_UU_ww * Aww_TF;

rhs.chi = advec.chi + two_over_gr_spacedim * vars.chi *
                      (vars.lapse * vars.K - divshift);

FOR(i, j)
{
rhs.h[i][j] = advec.h[i][j] - 2.0 * vars.lapse * A_TF[i][j] -
two_over_gr_spacedim * vars.h[i][j] * divshift;
FOR(k)
        {
                rhs.h[i][j] +=
                        vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
        }
}

rhs.hww = advec.hww - 2.0 * vars.lapse * Aww_TF -
          two_over_gr_spacedim * vars.hww * divshift_w;

Tensor<2, data_t> Adot_TF;
FOR(i, j)
{
Adot_TF[i][j] =
-covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
}
data_t Adot_TF_ww = -covd2lapse_ww + vars.chi * vars.lapse * ricci.LLww;

data_t trAdot = compute_trace(Adot_TF, h_UU) + nS * h_UU_ww * Adot_TF_ww;
// make_trace_free(Adot_TF, vars.h, h_UU);
FOR(i, j) { Adot_TF[i][j] -= one_over_gr_spacedim * vars.h[i][j] * trAdot; }
Adot_TF_ww -= one_over_gr_spacedim * vars.hww * trAdot;

FOR(i, j)
{
rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] +
vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) -
two_over_gr_spacedim * divshift);
FOR(k)
        {
                rhs.A[i][j] +=
                        A_TF[k][i] * d1.shift[k][j] + A_TF[k][j] * d1.shift[k][i];
        FOR(l)
        {
            rhs.A[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * A_TF[i][k] * A_TF[l][j];
        }
        }
}

rhs.Aww = advec.Aww + Adot_TF_ww +
          Aww_TF * (vars.lapse *
                    ((vars.K - 2 * vars.Theta) - 2 * h_UU_ww * Aww_TF) -
                    two_over_gr_spacedim * divshift_w);

#ifdef COVARIANTZ4
data_t kappa1_lapse = this->m_params.kappa1;
#else
data_t kappa1_lapse = this->m_params.kappa1 * vars.lapse;
#endif

rhs.Theta =
advec.Theta +
0.5 * vars.lapse *
(ricci.scalar - tr_A2 +
 ((GR_SPACEDIM - 1.0) * one_over_gr_spacedim) * vars.K * vars.K -
 2 * vars.Theta * vars.K) -
0.5 * vars.Theta * kappa1_lapse *
((GR_SPACEDIM + 1) + this->m_params.kappa2 * (GR_SPACEDIM - 1)) -
Z_dot_d1lapse;

rhs.Theta += -vars.lapse * this->m_cosmological_constant;

rhs.K =
advec.K +
vars.lapse * (ricci.scalar + vars.K * (vars.K - 2 * vars.Theta)) -
kappa1_lapse * GR_SPACEDIM * (1 + this->m_params.kappa2) * vars.Theta -
tr_covd2lapse;
rhs.K += -2 * vars.lapse * GR_SPACEDIM / ((double)GR_SPACEDIM - 1.) *
this->m_cosmological_constant;

Tensor<1, data_t> Gammadot;
FOR(i)
        {
                Gammadot[i] =
                        two_over_gr_spacedim *
                        (divshift * (chris_contracted[i] +
                                     2 * this->m_params.kappa3 * Z_over_chi[i]) -
                         2 * vars.lapse * vars.K * Z_over_chi[i]) -
                        2 * kappa1_lapse * Z_over_chi[i] +
                        nS *
                        (2. * vars.lapse * chris_ww[i] * A_UU_ww + h_UU_ww * reg_03[i]);
        FOR(j)
        {
            Gammadot[i] +=
                    2 * h_UU[i][j] *
                    (vars.lapse * d1.Theta[j] - vars.Theta * d1.lapse[j]) -
                    2 * A_UU[i][j] * d1.lapse[j] -
                    vars.lapse * (two_over_gr_spacedim * (GR_SPACEDIM - 1.0) *
                                  h_UU[i][j] * d1.K[j] +
                                  GR_SPACEDIM * A_UU[i][j] * d1.chi[j] / vars.chi) -
                    (chris_contracted[j] +
                     2 * this->m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j] +
                    (nS * (GR_SPACEDIM - 2.0) * one_over_gr_spacedim) * h_UU[i][j] *
                    reg_04[j];
            FOR(k)
            {
                Gammadot[i] +=
                        2 * vars.lapse * chris.ULL[i][j][k] * A_UU[j][k] +
                        h_UU[j][k] * d2.shift[i][j][k] +
                        ((GR_SPACEDIM - 2.0) * one_over_gr_spacedim) * h_UU[i][j] *
                        d2.shift[k][j][k];
            }
        }
        rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
        }

this->m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);
}

// Adds the matter source terms of the CCZ4 equations for a given 3+1 decomposition of the stress tensor (moved verbatim
// from CCZ4Cartoon::rhs_equation)
template <class gauge_t, class deriv_t>
template <class data_t, template <typename> class vars_t>
void CCZ4CartoonBase<gauge_t, deriv_t>::add_matter_source_terms(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                                                               const emtensorCartoon_t<data_t> &emtensor) const
{
using namespace TensorAlgebra;
const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
auto h_UU = compute_inverse_sym(vars.h);

Tensor<2, data_t> Sij_TF;
data_t Sww_TF;

FOR(i, j)
{
// Need to divide by chi here to compensate for factor of chi in
// emtensor.S
Sij_TF[i][j] = emtensor.Sij[i][j] - one_over_gr_spacedim * emtensor.S *
        vars.h[i][j] / vars.chi;
}
Sww_TF =
emtensor.Sww - one_over_gr_spacedim * emtensor.S * vars.hww / vars.chi;

FOR(i, j)
{
rhs.A[i][j] +=
-8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sij_TF[i][j];
}

rhs.Aww += -8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sww_TF;

rhs.K +=
4.0 * M_PI * m_G_Newton * vars.lapse * (emtensor.S - 3 * emtensor.rho);
rhs.Theta += -8.0 * M_PI * m_G_Newton * vars.lapse * emtensor.rho;

FOR(i)
        {
                data_t matter_term_Gamma = 0.0;
        FOR(j)
        {
            matter_term_Gamma += -16.0 * M_PI * m_G_Newton * vars.lapse *
                                 h_UU[i][j] * emtensor.Si[j];
        }

        rhs.Gamma[i] += matter_term_Gamma;
        rhs.B[i] += matter_term_Gamma;
        }
}

#endif /* CCZ4CARTOONBASE_IMPL_HPP_ */
