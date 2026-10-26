/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4CARTOONBASE_HPP_
#define CCZ4CARTOONBASE_HPP_

#include "CCZ4CartoonGeometry.hpp"
#include "CCZ4CartoonVars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

/// Matter-independent part of the CCZ4 cartoon right-hand side
/**
 * The vacuum CCZ4 equations with the cartoon terms and the gauge, and the coupling of a given 3+1 decomposition of the
 * stress tensor (emtensorCartoon_t) to the CCZ4 variables. CCZ4Cartoon (complex scalar field) and CCZ4CartoonDirac
 * (Dirac field) derive from this class and add their matter.
 **/
template <class gauge_t = MovingPunctureGauge, class deriv_t = FourthOrderDerivatives>
class CCZ4CartoonBase : public CCZ4RHS<gauge_t, deriv_t>
{
  public:
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;
    using params_t = typename CCZ4RHS<gauge_t, deriv_t>::params_t;

    CCZ4CartoonBase(params_t params, //!< The CCZ4 cartoon parameters
                    double dx,       //!< The grid spacing
                    double sigma,    //!< Kreiss-Oliger dissipation coefficient
                    double a_G_Newton,
                    int formulation = CCZ4::USE_CCZ4, //!< Switches between CCZ4, BSSN,...
                    double cosmological_constant = 0  //!< Value of the cosmological const.
    );

  protected:
    /// Vacuum CCZ4 rhs with the cartoon terms (and the gauge); the matter fields' own rhs are left untouched
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    void vacuum_rhs_equation(vars_t<data_t> &rhs, const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
                             const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
                             const double &cartoon_coord) const;

    /// Adds the matter terms (rho, S_i, S_ij, S_ww, S) to the rhs of A_ij, A_ww, K, Theta, Gamma^i and B^i
    template <class data_t, template <typename> class vars_t>
    void add_matter_source_terms(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                                 const emtensorCartoon_t<data_t> &emtensor) const;

    const double m_G_Newton;
};

#include "CCZ4CartoonBase.impl.hpp"

#endif /* CCZ4CARTOONBASE_HPP_ */
