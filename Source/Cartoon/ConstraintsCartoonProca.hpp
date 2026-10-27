/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CONSTRAINTSCARTOONPROCA_HPP_
#define CONSTRAINTSCARTOONPROCA_HPP_

#include "CCZ4CartoonProca.hpp"
#include "ConstraintsCartoonBase.hpp"

/// Hamiltonian and momentum constraints with a Proca field, plus the Proca diagnostics
/**
 * Stores Ham, Mom1, Mom2, Mom, rho, rho_ADM, Sx, Sy (as Constraints does for the scalar field) and
 *   N        = n chi^{-3/2} 2 pi y, the integrand of the particle number (sum it as c_N is summed for the boson stars:
 *              N = int n sqrt(gamma) d^3x),
 *   mod_Xphi = |X_phi|, which equals the 1D solver's central amplitude A_central = f(0)/alpha(0) at the centre of a star
 *              at rest (so max mod_Xphi tracks the central amplitude like max mod_phi for boson stars),
 *   n        = Im(E^i Xbar_i), the particle number density measured by the normal observers (n = -n_a j^a),
 *   Gauss    = |D_i E^i + mu^2 X_phi|, the modulus of the complex Gauss constraint.
 * The Proca stress tensor needs no lapse or shift, but the gauge variables are loaded with the rest.
 **/
class ConstraintsCartoonProca : public ConstraintsCartoonBase
{
  public:
    template <class data_t> using Vars = CCZ4CartoonProcaVars::VarsWithGauge<data_t>;
    template <class data_t> using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;

    ConstraintsCartoonProca(double dx, ProcaFieldCartoon::params_t a_proca_params, double a_G_Newton,
                            double cosmological_constant = 0)
        : ConstraintsCartoonBase(dx, a_G_Newton, cosmological_constant), m_proca(a_proca_params)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        Coordinates<data_t> coords{current_cell, m_dx};

        CartoonConstraints_t<data_t> out = vacuum_constraint_equations(vars, d1, d2, coords.y);

        ProcaSources_t<data_t> src;
        m_proca.compute_sources(src, vars, d1, coords.y);
        add_matter_constraint_terms(out, vars, proca_emtensor_cartoon(src, vars));

        current_cell.store_vars(out.Ham, c_Ham);
        current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom2>());
        current_cell.store_vars(out.Mom_abs, c_Mom);
        current_cell.store_vars(out.rho, c_rho);
        current_cell.store_vars(out.rho_ADM, c_rho_ADM);
        current_cell.store_vars(out.Si[0], c_Sx);
        current_cell.store_vars(out.Si[1], c_Sy);

        const data_t N = src.n * pow(vars.chi, -1.5) * 2. * M_PI * coords.y;
        current_cell.store_vars(N, c_N);
        current_cell.store_vars(sqrt(vars.Xphi_Re * vars.Xphi_Re + vars.Xphi_Im * vars.Xphi_Im), c_mod_Xphi);
        current_cell.store_vars(src.n, c_n);
        current_cell.store_vars(sqrt(src.C_Re * src.C_Re + src.C_Im * src.C_Im), c_Gauss);
    }

  protected:
    ProcaFieldCartoon m_proca;
};

#endif /* CONSTRAINTSCARTOONPROCA_HPP_ */
