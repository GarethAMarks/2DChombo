/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CONSTRAINTSCARTOONBASE_HPP_
#define CONSTRAINTSCARTOONBASE_HPP_

#include "CCZ4CartoonGeometry.hpp"
#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FArrayBox.H"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

//! Output of the cartoon constraint calculation
template <class data_t> struct CartoonConstraints_t
{
    data_t Ham;
    Tensor<1, data_t> Mom;
    data_t Mom_abs;

    data_t rho;
    data_t rho_ADM;
    Tensor<1, data_t> Si;
    Tensor<2, data_t> Sij;
    data_t Sww;
    data_t S;
    Tensor<2, data_t> Sij_TF;
    data_t Sww_TF;
};

/// Matter-independent part of the cartoon Hamiltonian and momentum constraints
/**
 * Constraints (complex scalar field) and ConstraintsCartoonDirac (Dirac field) derive from this class and add their
 * matter through add_matter_constraint_terms.
 **/
class ConstraintsCartoonBase
{
  public:
    ConstraintsCartoonBase(double dx, double a_G_Newton, double cosmological_constant = 0);

  protected:
    const FourthOrderDerivatives m_deriv;
    double m_cosmological_constant;
    double m_dx;
    double m_G_Newton;

    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    CartoonConstraints_t<data_t> vacuum_constraint_equations(const vars_t<data_t> &vars,
                                                             const vars_t<Tensor<1, data_t>> &d1,
                                                             const diff2_vars_t<Tensor<2, data_t>> &d2,
                                                             const double &cartoon_coord) const;

    //! adds -16 pi G rho to Ham and -8 pi G S_i to Mom, then fills Mom_abs, rho, rho_ADM and Si
    template <class data_t, template <typename> class vars_t>
    void add_matter_constraint_terms(CartoonConstraints_t<data_t> &out, const vars_t<data_t> &vars,
                                     const emtensorCartoon_t<data_t> &emtensor) const;
};

#include "ConstraintsCartoonBase.impl.hpp"

#endif /* CONSTRAINTSCARTOONBASE_HPP_ */
