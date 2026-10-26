/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTSCARTOONBASE_HPP_)
#error "This file should only be included through ConstraintsCartoonBase.hpp"
#endif

#ifndef CONSTRAINTSCARTOONBASE_IMPL_HPP_
#define CONSTRAINTSCARTOONBASE_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

inline ConstraintsCartoonBase::ConstraintsCartoonBase(double dx, double a_G_Newton, double cosmological_constant)
    : m_deriv(dx), m_cosmological_constant(cosmological_constant), m_dx(dx), m_G_Newton(a_G_Newton)
{
}

// Vacuum Hamiltonian and momentum constraints with the cartoon terms (moved verbatim from Constraints::constraint_equations)
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
CartoonConstraints_t<data_t> ConstraintsCartoonBase::vacuum_constraint_equations(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
    const double &cartoon_coord) const
{
    using namespace TensorAlgebra;

    CartoonConstraints_t<data_t> out;

    const int dI = CH_SPACEDIM - 1;
    const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
    const double one_over_cartoon_coord = 1. / cartoon_coord;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    auto h_UU_ww = 1. / vars.hww;
    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    Tensor<1, data_t> chris_ww;
    FOR(i)
    {
        chris_ww[i] =
            one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
    }

    Tensor<1, data_t> reg_07;
    FOR(i)
    {
        reg_07[i] =
            one_over_cartoon_coord * (vars.A[i][dI] - delta(i, dI) * vars.Aww);
    }

    auto ricci = CCZ4CartoonGeometry::compute_ricci(vars, d1, d2, h_UU, h_UU_ww,
                                                    chris, cartoon_coord);

    auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
    data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU) +
                   nS * vars.Aww * vars.Aww * h_UU_ww * h_UU_ww;

    out.Ham = ricci.scalar +
              (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;
    out.Ham -= 2 * m_cosmological_constant;

    Tensor<2, data_t> covd_A[CH_SPACEDIM];
    FOR(i, j, k)
    {
        covd_A[i][j][k] = d1.A[j][k][i];
        FOR(l)
        {
            covd_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                               chris.ULL[l][i][k] * vars.A[l][j];
        }
    }

    FOR(i)
    {
        out.Mom[i] =
            -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM +
            nS * h_UU_ww * (reg_07[i] - 0.5 * h_UU_ww * vars.Aww * d1.hww[i]);
    }
    FOR(i, j)
    {
        out.Mom[i] -= nS * h_UU_ww * chris_ww[j] * vars.A[i][j];
        FOR(k)
        {
            out.Mom[i] += h_UU[j][k] * (covd_A[k][j][i] -
                                        GR_SPACEDIM * vars.A[i][j] * d1.chi[k] /
                                            (2 * chi_regularised));
        }
    }

    return out;
}

// Matter terms of the constraints and the stored source diagnostics (moved verbatim from Constraints::constraint_equations)
template <class data_t, template <typename> class vars_t>
void ConstraintsCartoonBase::add_matter_constraint_terms(CartoonConstraints_t<data_t> &out, const vars_t<data_t> &vars,
                                                         const emtensorCartoon_t<data_t> &emtensor) const
{
    out.Ham += -16.0 * M_PI * m_G_Newton * emtensor.rho;
    FOR(i)
    {
        out.Mom[i] += -8.0 * M_PI * m_G_Newton * emtensor.Si[i];
    }

    data_t Mom_abs = 0.;

    FOR(i) { Mom_abs += out.Mom[i] * out.Mom[i]; }
    Mom_abs = sqrt(Mom_abs);
    out.Mom_abs = Mom_abs;

    out.rho = emtensor.rho;
    data_t det_gamma = TensorAlgebra::compute_determinant(vars.h);
    det_gamma *= vars.hww;
    out.rho_ADM = out.rho * det_gamma;
    FOR(i)
    {
        out.Si[i] = emtensor.Si[i];
    }
}

#endif /* CONSTRAINTSCARTOONBASE_IMPL_HPP_ */
