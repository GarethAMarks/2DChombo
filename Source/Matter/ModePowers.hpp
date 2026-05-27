/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODEPOWERS_HPP_
#define MODEPOWERS_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the real and imaginary contributions to the mode coefficients
//! of |phi| exp(i m psi) for m = 0..12 on the grid.
template <class deriv_t = FourthOrderDerivatives>
class ModePowers
{
  protected:
    deriv_t m_deriv;

    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;

  public:
    ModePowers(const double a_dx) : m_deriv(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

        data_t mod_phi = sqrt(vars.phi * vars.phi + vars.phi_Im * vars.phi_Im);
        data_t psi = atan2(coords.y, coords.x);

        for (int m = 0; m <= 12; ++m)
        {
            data_t angle = m * psi;
            current_cell.store_vars(mod_phi * cos(angle), c_mode_power_re_0 + m);
            current_cell.store_vars(mod_phi * sin(angle), c_mode_power_im_0 + m);
        }
    }
};

#endif /* MODEPOWERS_HPP_ */
