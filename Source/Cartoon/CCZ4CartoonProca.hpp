/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4CARTOONPROCA_HPP_
#define CCZ4CARTOONPROCA_HPP_

#include "CCZ4CartoonBase.hpp"
#include "ProcaFieldCartoon.hpp"

/// CCZ4 cartoon variables plus the Proca field (X_x, X_y, E^x, E^y, X_phi, Z, real and imaginary parts)
namespace CCZ4CartoonProcaVars
{
template <class data_t>
struct VarsWithGauge : public CCZ4CartoonVars::VarsWithGauge<data_t>, public ProcaCartoonVars<data_t>
{
    template <typename mapping_function_t> void enum_mapping(mapping_function_t mapping_function)
    {
        CCZ4CartoonVars::VarsWithGauge<data_t>::enum_mapping(mapping_function);
        ProcaCartoonVars<data_t>::enum_mapping(mapping_function);
    }
};

template <class data_t>
struct Diff2VarsWithGauge : public CCZ4CartoonVars::Diff2VarsWithGauge<data_t>, public ProcaCartoonDiff2Vars<data_t>
{
    template <typename mapping_function_t> void enum_mapping(mapping_function_t mapping_function)
    {
        CCZ4CartoonVars::Diff2VarsWithGauge<data_t>::enum_mapping(mapping_function);
        ProcaCartoonDiff2Vars<data_t>::enum_mapping(mapping_function);
    }
};
} // namespace CCZ4CartoonProcaVars

/// Converts the Proca sources on the cartoon plane to the emtensorCartoon_t used by the CCZ4 cartoon equations
template <class data_t, template <typename> class vars_t>
emtensorCartoon_t<data_t> proca_emtensor_cartoon(const ProcaSources_t<data_t> &src, const vars_t<data_t> &vars)
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

/// Compute class for the CCZ4 cartoon right-hand side with a complex Proca field (see ProcaFieldCartoon.hpp)
/**
 * The vacuum part is CCZ4CartoonBase (shared with the scalar field version CCZ4Cartoon and the Dirac version); the
 * Proca field adds its stress tensor to the CCZ4 equations and is evolved with the Zilhao-Witek-Cardoso equations.
 * Kreiss-Oliger dissipation is applied to all variables, including the Proca field (the first-order Proca system
 * needs it).
 **/
template <class gauge_t = MovingPunctureGauge, class deriv_t = FourthOrderDerivatives>
class CCZ4CartoonProca : public CCZ4CartoonBase<gauge_t, deriv_t>
{
  public:
    using params_t = typename CCZ4CartoonBase<gauge_t, deriv_t>::params_t;
    using CCZ4 = typename CCZ4CartoonBase<gauge_t, deriv_t>::CCZ4;

    template <class data_t> using Vars = CCZ4CartoonProcaVars::VarsWithGauge<data_t>;
    template <class data_t> using Diff2Vars = CCZ4CartoonProcaVars::Diff2VarsWithGauge<data_t>;

    CCZ4CartoonProca(params_t params, double dx, double sigma, ProcaFieldCartoon::params_t a_proca_params,
                     double a_G_Newton, int formulation = CCZ4::USE_CCZ4, double cosmological_constant = 0)
        : CCZ4CartoonBase<gauge_t, deriv_t>(params, dx, sigma, a_G_Newton, formulation, cosmological_constant),
          m_proca(a_proca_params)
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

        ProcaCartoonVars<data_t> rhs_p;
        ProcaSources_t<data_t> src;
        m_proca.compute(rhs_p, src, vars, d1, d2, coords.y);
        for (int i = 0; i < 2; ++i)
        {
            rhs.X_Re[i] = rhs_p.X_Re[i] + advec.X_Re[i];
            rhs.X_Im[i] = rhs_p.X_Im[i] + advec.X_Im[i];
            rhs.E_Re[i] = rhs_p.E_Re[i] + advec.E_Re[i];
            rhs.E_Im[i] = rhs_p.E_Im[i] + advec.E_Im[i];
        }
        rhs.Xphi_Re = rhs_p.Xphi_Re + advec.Xphi_Re;
        rhs.Xphi_Im = rhs_p.Xphi_Im + advec.Xphi_Im;
        rhs.Z_Re = rhs_p.Z_Re + advec.Z_Re;
        rhs.Z_Im = rhs_p.Z_Im + advec.Z_Im;

        this->add_matter_source_terms(rhs, vars, proca_emtensor_cartoon(src, vars));

        this->m_deriv.add_dissipation(rhs, current_cell, this->m_sigma);

        current_cell.store_vars(rhs);
    }

  protected:
    ProcaFieldCartoon m_proca;
};

#endif /* CCZ4CARTOONPROCA_HPP_ */
