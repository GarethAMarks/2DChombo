/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CONSTRAINTSCARTOON_HPP_)
#error "This file should only be included through ConstraintsCartoon.hpp"
#endif

#ifndef CONSTRAINTSCARTOON_IMPL_HPP_
#define CONSTRAINTSCARTOON_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

template <class potential_t>
inline Constraints<potential_t>::Constraints(
    double dx, potential_t a_potential, double a_G_Newton,
    double cosmological_constant /*defaulted*/)
    : ConstraintsCartoonBase(dx, a_G_Newton, cosmological_constant),
      m_potential(a_potential)
{
}

template <class potential_t>
template <class data_t>
void Constraints<potential_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    Coordinates<data_t> coords{current_cell, m_dx};

    // Function below calls function to compute matter terms, too
    constraints_t<data_t> out = constraint_equations(vars, d1, d2, coords.y);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Ham, c_Ham);
#if CH_SPACEDIM == 3
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom3>());
#elif CH_SPACEDIM == 2
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom2>());
#else
#ifdef CH_SPACEDIM
#error current_cell.store_vars() has not got your dimension combination implemented.
#endif
#endif

    current_cell.store_vars(out.Mom_abs, c_Mom);
    current_cell.store_vars(out.rho, c_rho);
    current_cell.store_vars(out.rho_ADM, c_rho_ADM);
    current_cell.store_vars(out.Si[0], c_Sx);
    current_cell.store_vars(out.Si[1], c_Sy);
}

template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
typename Constraints<potential_t>:: template constraints_t<data_t>
Constraints<potential_t>::constraint_equations(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const double &cartoon_coord) const
{
    // vacuum part (ConstraintsCartoonBase), then the complex scalar field
    constraints_t<data_t> out =
        this->vacuum_constraint_equations(vars, d1, d2, cartoon_coord);

    using namespace TensorAlgebra;
    const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    auto h_UU_ww = 1. / vars.hww;
    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Matter contributions
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi, vars);
    
    emtensorCartoon_t<data_t> emtensor = compute_SF_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, V_of_phi);
    this->add_matter_constraint_terms(out, vars, emtensor);

    return out;
}

#endif /* CONSTRAINTSCARTOON_IMPL_HPP_ */
