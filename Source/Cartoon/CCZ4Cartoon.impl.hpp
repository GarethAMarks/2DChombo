/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(CCZ4CARTOON_HPP_)
#error "This file should only be included through CCZ4Cartoon.hpp"
#endif

#ifndef CCZ4CARTOON_IMPL_HPP_
#define CCZ4CARTOON_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"
#include <iostream>
#include <stdio.h>

template <class gauge_t, class deriv_t, class potential_t>
inline CCZ4Cartoon<gauge_t, deriv_t, potential_t>::CCZ4Cartoon(
        params_t params, double dx, double sigma,
        potential_t a_potential, double a_G_Newton, int formulation,
        double cosmological_constant)
        : CCZ4CartoonBase<gauge_t, deriv_t>(params, dx, sigma, a_G_Newton,
                                            formulation, cosmological_constant),
          m_potential(a_potential)
{
}

template <class data_t, template <typename> class vars_t>
emtensorCartoon_t<data_t>
compute_SF_EM_tensor(const vars_t<data_t> &vars,
                     const vars_t<Tensor<1, data_t>> &d1,
const Tensor<2, data_t> &h_UU, const data_t &h_UU_ww,
const chris_t<data_t> &chris, const int &nS,
const data_t V_of_phi)
{
emtensorCartoon_t<data_t> out;

// Useful quantity Vt
// 1st term of rho
data_t Vt = vars.Pi * vars.Pi + vars.Pi_Im * vars.Pi_Im;
//std::cout << "vars.Pi"<<vars.Pi<<std::endl;
//std::cout << "vars.Pi_Im"<<vars.Pi_Im<<std::endl;
//std::cout << "vars.phi"<<vars.phi<<std::endl;
//std::cout << "vars.phi_Im"<<vars.phi_Im<<std::endl;
// 2nd term of rho
//std::cout << "Vt"<<Vt<<std::endl;


data_t Vt1 = 0;
FOR(i, j) { Vt1 += vars.chi * h_UU[i][j] * (d1.phi[i] * d1.phi[j] + d1.phi_Im[i] * d1.phi_Im[j]); }
//std::cout << "Vt1"<<Vt1<<std::endl;



// Full rho (added potential already)
out.rho = 0.5 * (Vt + Vt1) + 0.5 * V_of_phi;
//std::cout << "out.rho"<<out.rho<<std::endl;
// S_i
FOR(i) { out.Si[i] = vars.Pi * d1.phi[i] + vars.Pi_Im * d1.phi_Im[i];

//std::cout << "out.Si"<<i<<out.Si[i]<<std::endl;
}

// S_ij
FOR(i, j)
{
out.Sij[i][j] =
-0.5 * vars.h[i][j] * (-Vt + Vt1 + V_of_phi) / vars.chi + d1.phi[i] * d1.phi[j] + d1.phi_Im[i] * d1.phi_Im[j];

//std::cout << "out.Sij"<<i<<j<<out.Sij[i][j]<<std::endl;
}



// Sww
out.Sww = -0.5 * vars.hww * (-Vt + Vt1 + V_of_phi) / vars.chi;

// S
out.S = TensorAlgebra::compute_trace(out.Sij, h_UU);
out.S += nS * h_UU_ww * out.Sww;
out.S *= vars.chi;

//std::cout << "out.S"<<out.S<<std::endl;
// Potential terms, I hate this way!!!
//out.rho += 0.5 * V_of_phi;
//out.S += -1.5 * V_of_phi;

// NEED TO CHECK THESE BELOW!!!
//FOR(i, j)
//{
//out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi;
//}
//out.Sww += -0.5 * vars.hww * V_of_phi / vars.chi;

return out;
}

template <class gauge_t, class deriv_t, class potential_t>
template <class data_t>
void CCZ4Cartoon<gauge_t, deriv_t, potential_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

    Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

    Vars<data_t> rhs;
    this->rhs_equation(rhs, vars, d1, d2, advec, coords.y);

    this->m_deriv.add_dissipation(rhs, current_cell, this->m_sigma);

    current_cell.store_vars(rhs); // Write the rhs into the output FArrayBox
}

template <class gauge_t, class deriv_t, class potential_t>
template <class data_t, template <typename> class vars_t,
        template <typename> class diff2_vars_t>
void CCZ4Cartoon<gauge_t, deriv_t, potential_t>::rhs_equation(
        vars_t<data_t> &rhs, const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1,
const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
const double &cartoon_coord) const
{
// vacuum CCZ4 cartoon terms and gauge (CCZ4CartoonBase), then the complex scalar field
this->vacuum_rhs_equation(rhs, vars, d1, d2, advec, cartoon_coord);

using namespace TensorAlgebra;
const int dI = CH_SPACEDIM - 1;
const int nS =
        GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
auto h_UU = compute_inverse_sym(vars.h);
auto h_UU_ww = 1. / vars.hww;
auto chris = compute_christoffel(d1.h, h_UU);
const double one_over_cartoon_coord = 1. / cartoon_coord;

// Matter contributions

// Need to write a function that returns potential at some point.
data_t V_of_phi = 0.0;  // note that! here is V_of_modulus_phi_squared actually. but i am lazy to modify
data_t dVdmodulus_phi_squared = 0.0;

m_potential.compute_potential(V_of_phi, dVdmodulus_phi_squared, vars);

emtensorCartoon_t<data_t> emtensor =
        compute_SF_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, V_of_phi);

this->add_matter_source_terms(rhs, vars, emtensor);

// For matter fields
rhs.phi = -vars.lapse * vars.Pi + advec.phi;
rhs.phi_Im = -vars.lapse * vars.Pi_Im + advec.phi_Im;
//std::cout << "rhs.phi"<<rhs.phi<<std::endl;
//std::cout << "rhs.phi_Im"<<rhs.phi_Im<<std::endl;

// The Pi_Re part
// 7th term + 1st term + 6th term + potential term
rhs.Pi = vars.lapse * vars.K * vars.Pi
         + advec.Pi
         -vars.lapse * vars.chi * h_UU_ww * d1.phi[dI] * one_over_cartoon_coord
         +vars.lapse * dVdmodulus_phi_squared * vars.phi;

// 2nd term + 4th term +  5th term
FOR(i, j)
{
// includes non conformal parts of chris not included in chris_ULL
rhs.Pi += h_UU[i][j] * (-vars.chi * d1.lapse[j] * d1.phi[i]
+0.5 * d1.chi[i] * vars.lapse * d1.phi[j]
-vars.chi * vars.lapse * d2.phi[i][j]);

}


// the 3rd term
FOR(i) { rhs.Pi += vars.chi * vars.lapse * vars.Gamma[i] * d1.phi[i]; }
// the 8th term
//rhs.Pi += vars.lapse * dVdphi;
//std::cout << "rhs.Pi"<<rhs.Pi<<std::endl;


//////////////////////////////////////im part


// The Pi_Im part
// 7th term + 1st term + 6th term + potential term
rhs.Pi_Im = vars.lapse * vars.K * vars.Pi_Im
            + advec.Pi_Im
            - vars.lapse * vars.chi * h_UU_ww * d1.phi_Im[dI] * one_over_cartoon_coord
            + vars.lapse * dVdmodulus_phi_squared * vars.phi_Im;

// 2nd term + 4th term +  5th term
FOR(i, j)
{
// includes non conformal parts of chris not included in chris_ULL
rhs.Pi_Im += h_UU[i][j] * (-vars.chi * d1.lapse[j] * d1.phi_Im[i]
             + 0.5 * d1.chi[i] * vars.lapse * d1.phi_Im[j]
             - vars.chi * vars.lapse * d2.phi_Im[i][j]);

}


// the 3rd term
FOR(i) { rhs.Pi_Im += vars.chi * vars.lapse * vars.Gamma[i] * d1.phi_Im[i]; }

//std::cout << "rhs.Pi_Im"<<rhs.Pi_Im<<std::endl;
// the 8th term
//rhs.Pi_Im += vars.lapse * dVdphi_Im;



}

#endif /* CCZ4CARTOON_IMPL_HPP_ */
