/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CONSTRAINTSCARTOONDIRAC_HPP_
#define CONSTRAINTSCARTOONDIRAC_HPP_

#include "CCZ4CartoonDirac.hpp"
#include "ConstraintsCartoonBase.hpp"

/// Hamiltonian and momentum constraints with Dirac matter, plus the Dirac diagnostics
/**
 * Stores Ham, Mom1, Mom2, Mom, rho, rho_ADM, Sx, Sy (as Constraints does for the scalar field) and
 *   N       = n chi^{-3/2} 2 pi y,  the integrand of the particle number of both spinors (sum it as c_N is summed for the
 *             boson stars: N = int n sqrt(gamma) d^3x),
 *   mod_psi = sqrt(4 pi psi^dag psi), which equals the 1D solver's central amplitude A_central = g(0) at the centre of
 *             the unperturbed star (so max mod_psi tracks the central amplitude like max mod_phi for boson stars),
 *   n       = the particle number density of both spinors measured by the normal observers (n = -n_mu J^mu; the
 *             particle number is int n sqrt(gamma) d^3x),
 *   n_psi   = psi^dag psi, the number density of each spinor field: the evolved J_x = +1/2 spinor and its mirror image
 *             (the J_x = -1/2 spinor) have the same density on the cartoon plane, so n = 2 n_psi.
 * The energy density needs n^a nabla_a psi, which is eliminated with the Dirac equation, so the lapse and shift (and
 * K_ij) are loaded too.
 **/
class ConstraintsCartoonDirac : public ConstraintsCartoonBase
{
  public:
    template <class data_t> using Vars = CCZ4CartoonDiracVars::VarsWithGauge<data_t>;
    template <class data_t> using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;

    ConstraintsCartoonDirac(double dx, DiracFieldCartoon::params_t a_dirac_params, double a_G_Newton,
                            double cosmological_constant = 0)
        : ConstraintsCartoonBase(dx, a_G_Newton, cosmological_constant), m_dirac(a_dirac_params)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
        Coordinates<data_t> coords{current_cell, m_dx};

        CartoonConstraints_t<data_t> out = vacuum_constraint_equations(vars, d1, d2, coords.y);

        DiracSpinorVars<data_t> rhs_psi;
        DiracSources_t<data_t> src;
        m_dirac.compute(rhs_psi, src, vars, d1, coords.y);
        add_matter_constraint_terms(out, vars, dirac_emtensor_cartoon(src, vars));

        current_cell.store_vars(out.Ham, c_Ham);
        current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom2>());
        current_cell.store_vars(out.Mom_abs, c_Mom);
        current_cell.store_vars(out.rho, c_rho);
        current_cell.store_vars(out.rho_ADM, c_rho_ADM);
        current_cell.store_vars(out.Si[0], c_Sx);
        current_cell.store_vars(out.Si[1], c_Sy);

        const data_t N = src.n * pow(vars.chi, -1.5) * 2. * M_PI * coords.y;
        current_cell.store_vars(N, c_N);
        current_cell.store_vars(sqrt(2. * M_PI * src.n), c_mod_psi);
        current_cell.store_vars(src.n, c_n);
        current_cell.store_vars(0.5 * src.n, c_n_psi);
    }

  protected:
    DiracFieldCartoon m_dirac;
};

#endif /* CONSTRAINTSCARTOONDIRAC_HPP_ */
