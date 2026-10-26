/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4CARTOONDIRAC_HPP_
#define CCZ4CARTOONDIRAC_HPP_

#include "CCZ4CartoonBase.hpp"
#include "DiracFieldCartoon.hpp"

/// CCZ4 cartoon variables plus the Dirac spinor (psi0 ... psi3, real and imaginary parts)
namespace CCZ4CartoonDiracVars
{
template <class data_t>
struct VarsWithGauge : public CCZ4CartoonVars::VarsWithGauge<data_t>, public DiracSpinorVars<data_t>
{
    template <typename mapping_function_t> void enum_mapping(mapping_function_t mapping_function)
    {
        CCZ4CartoonVars::VarsWithGauge<data_t>::enum_mapping(mapping_function);
        DiracSpinorVars<data_t>::enum_mapping(mapping_function);
    }
};
} // namespace CCZ4CartoonDiracVars

/// Converts the Dirac sources on the cartoon plane to the emtensorCartoon_t used by the CCZ4 cartoon equations
template <class data_t, template <typename> class vars_t>
emtensorCartoon_t<data_t> dirac_emtensor_cartoon(const DiracSources_t<data_t> &src, const vars_t<data_t> &vars)
{
    emtensorCartoon_t<data_t> out;
    const int nS = GR_SPACEDIM - CH_SPACEDIM;
    out.rho = src.rho;
    for (int i = 0; i < 2; ++i)
    {
        out.Si[i] = src.Si[i];
        for (int j = 0; j < 2; ++j)
            out.Sij[i][j] = src.Sij[i][j];
    }
    out.Sww = src.Sww;
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    out.S = vars.chi * (TensorAlgebra::compute_trace(out.Sij, h_UU) + nS * out.Sww / vars.hww);
    return out;
}

/// Compute class for the CCZ4 cartoon right-hand side with a Dirac field (two spinors, see DiracFieldCartoon.hpp)
/**
 * The vacuum part is CCZ4CartoonBase (shared with the complex scalar field version CCZ4Cartoon); the Dirac field adds
 * its stress tensor to the CCZ4 equations and evolves psi. Kreiss-Oliger dissipation is applied to all variables,
 * including the spinor (as in the 1D code, the first-order Dirac system needs it).
 **/
template <class gauge_t = MovingPunctureGauge, class deriv_t = FourthOrderDerivatives>
class CCZ4CartoonDirac : public CCZ4CartoonBase<gauge_t, deriv_t>
{
  public:
    using params_t = typename CCZ4CartoonBase<gauge_t, deriv_t>::params_t;
    using CCZ4 = typename CCZ4CartoonBase<gauge_t, deriv_t>::CCZ4;

    template <class data_t> using Vars = CCZ4CartoonDiracVars::VarsWithGauge<data_t>;
    template <class data_t> using Diff2Vars = CCZ4CartoonVars::Diff2VarsWithGauge<data_t>;

    CCZ4CartoonDirac(params_t params, double dx, double sigma, DiracFieldCartoon::params_t a_dirac_params,
                     double a_G_Newton, int formulation = CCZ4::USE_CCZ4, double cosmological_constant = 0)
        : CCZ4CartoonBase<gauge_t, deriv_t>(params, dx, sigma, a_G_Newton, formulation, cosmological_constant),
          m_dirac(a_dirac_params)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
        const auto advec = this->m_deriv.template advection<Vars>(current_cell, vars.shift);

        Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

        Vars<data_t> rhs;
        this->vacuum_rhs_equation(rhs, vars, d1, d2, advec, coords.y);

        DiracSpinorVars<data_t> rhs_psi;
        DiracSources_t<data_t> src;
        m_dirac.compute(rhs_psi, src, vars, d1, coords.y);
        for (int c = 0; c < 4; ++c)
        {
            rhs.psi_Re[c] = rhs_psi.psi_Re[c] + advec.psi_Re[c];
            rhs.psi_Im[c] = rhs_psi.psi_Im[c] + advec.psi_Im[c];
        }

        this->add_matter_source_terms(rhs, vars, dirac_emtensor_cartoon(src, vars));

        this->m_deriv.add_dissipation(rhs, current_cell, this->m_sigma);

        current_cell.store_vars(rhs);
    }

  protected:
    DiracFieldCartoon m_dirac;
};

#endif /* CCZ4CARTOONDIRAC_HPP_ */
