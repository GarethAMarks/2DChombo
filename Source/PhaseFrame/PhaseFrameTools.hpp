/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHASEFRAMETOOLS_HPP_
#define PHASEFRAMETOOLS_HPP_

// Pointwise compute classes for the U(1) phase frame (use with BoxLoops):
//  PhaseFrameRotate       - rotate the matter slots by exp(+-i Theta)
//  PhaseFrameInitialOmega - initial frame frequency Omega(x, t=0)
//  PhaseFrameSetTestTheta - pure-gauge test profile Theta_0(x) (+ rotation)
//  PhaseFrameCopyPhysical - copy phi/phi_Im of a physical-frame state into
//                           the diagnostic variables phi_phys/phi_Im_phys
//
// Theta is an axisymmetric scalar (a rotation in field space), so it commutes
// with the cartoon reduction: nothing here depends on the reduction except
// that Theta(x, y) must be even in the cartoon coordinate y.

#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "PhaseFrameParams.hpp"
#include "PhaseFrameVars.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"
#include <cmath>

/// Rotates phi, Pi by exp(sign * i Theta) in place (pointwise).
/// sign = +1: frame (psi, P) -> physical (phi, Pi); sign = -1: the inverse.
class PhaseFrameRotate
{
  public:
    PhaseFrameRotate(const int a_sign = +1) : m_sign(a_sign) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const data_t theta = current_cell.load_vars(c_frame_theta);
        const auto in =
            current_cell.template load_vars<PhaseFrame::MatterVars>();
        const data_t c = cos(theta);
        const data_t s = static_cast<double>(m_sign) * sin(theta);

        PhaseFrame::MatterVars<data_t> out;
        PhaseFrame::complex_rotate(in.phi, in.phi_Im, c, s, out.phi,
                                   out.phi_Im);
        PhaseFrame::complex_rotate(in.Pi, in.Pi_Im, c, s, out.Pi, out.Pi_Im);
        current_cell.store_vars(out);
    }

  protected:
    const int m_sign;
};

/// Initial frame frequency. With omega_from_ID, the eps-weighted local
/// frequency of the (physical) initial data,
///   Omega = Im(conj(phi) d_t phi) / (|phi|^2 + eps^2),
///   d_t phi = -lapse Pi + shift^i d_i phi,
/// which is ~ the (boosted) frequency of each star inside it and -> 0 (no
/// rotation, i.e. the original variables) far from the stars. Otherwise the
/// uniform value omega0. Writes only frame_omega and frame_theta (= 0), so it
/// can be looped in place over the freshly set initial data.
class PhaseFrameInitialOmega
{
  public:
    PhaseFrameInitialOmega(const PhaseFrame_params_t &a_params,
                           const double a_dx)
        : m_params(a_params), m_deriv(a_dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        data_t omega = m_params.omega0;
        if (m_params.omega_from_ID)
        {
            const auto vars = current_cell.template load_vars<
                CCZ4CartoonVars::VarsWithGauge>();
            const auto advec =
                m_deriv.template advection<PhaseFrame::MatterVars>(
                    current_cell, vars.shift);
            const data_t dphi = -vars.lapse * vars.Pi + advec.phi;
            const data_t dphi_Im = -vars.lapse * vars.Pi_Im + advec.phi_Im;
            const data_t mod_phi2 =
                vars.phi * vars.phi + vars.phi_Im * vars.phi_Im;
            omega = (vars.phi * dphi_Im - vars.phi_Im * dphi) /
                    (mod_phi2 + m_params.eps * m_params.eps);
        }
        const data_t zero = 0.;
        current_cell.store_vars(zero, c_frame_theta);
        current_cell.store_vars(omega, c_frame_omega);
    }

  protected:
    const PhaseFrame_params_t m_params;
    const FourthOrderDerivatives m_deriv;
};

/// Pure-gauge test: sets Theta = a r^2/(r^2 + w^2), r^2 = (x - x_c)^2 + y^2,
/// centred on the symmetry axis so that Theta stays axisymmetric, and rotates
/// the (physical) matter slots into that frame, psi = exp(-i Theta) phi.
/// Pointwise, so it may be looped including ghost cells.
class PhaseFrameSetTestTheta
{
  public:
    PhaseFrameSetTestTheta(const PhaseFrame_params_t &a_params,
                           const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx);
        const data_t xr = coords.x - m_params.test_theta0_centre_x;
        const data_t r2 = xr * xr + coords.y * coords.y;
        const double w2 = m_params.test_theta0_width *
                          m_params.test_theta0_width;
        const data_t theta = m_params.test_theta0 * r2 / (r2 + w2);

        const auto in =
            current_cell.template load_vars<PhaseFrame::MatterVars>();
        const data_t c = cos(theta);
        const data_t s = -sin(theta);
        PhaseFrame::MatterVars<data_t> out;
        PhaseFrame::complex_rotate(in.phi, in.phi_Im, c, s, out.phi,
                                   out.phi_Im);
        PhaseFrame::complex_rotate(in.Pi, in.Pi_Im, c, s, out.Pi, out.Pi_Im);
        current_cell.store_vars(out);
        current_cell.store_vars(theta, c_frame_theta);
    }

  protected:
    const PhaseFrame_params_t m_params;
    const double m_dx;
};

/// Copies phi, phi_Im of a state that is already in the physical frame into
/// the diagnostic variables phi_phys, phi_Im_phys (for plot files).
class PhaseFrameCopyPhysical
{
  public:
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const data_t phi = current_cell.load_vars(c_phi);
        const data_t phi_Im = current_cell.load_vars(c_phi_Im);
        current_cell.store_vars(phi, c_phi_phys);
        current_cell.store_vars(phi_Im, c_phi_Im_phys);
    }
};

#endif /* PHASEFRAMETOOLS_HPP_ */
