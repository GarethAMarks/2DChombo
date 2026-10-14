/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHASEFRAMEPARAMS_HPP_
#define PHASEFRAMEPARAMS_HPP_

/// Parameters of the U(1) phase frame for the complex scalar field.
/**
 * The evolved matter variables are psi = exp(-i Theta) phi and
 * P = exp(-i Theta) Pi, with d_t Theta = Omega. See CCZ4CartoonPhaseFrame.hpp
 * for the equations and params.txt (ScalarField_2D) for a description of each
 * parameter.
 */
struct PhaseFrame_params_t
{
    bool use_phase_frame = false; //!< evolve psi, P, Theta, Omega

    /// true: link variables, i.e. every stencil acts on the physical field
    /// exp(i Theta_k) psi_k rotated into the frame of the central cell
    /// (semi-discretisation identical to the physical-field one, KO included).
    /// false: continuum connection D_i psi = d_i psi + i (d_i Theta) psi, with
    /// KO acting on psi directly.
    bool links = true;

    /// true: plain local driver for Omega; false: Omega frozen at its initial
    /// profile
    bool local_driver = true;

    /// initial Omega: true = local frequency of the initial data,
    /// Im(conj(phi) d_t phi)/(|phi|^2 + eps^2); false = uniform omega0
    bool omega_from_ID = true;
    double omega0 = 0.0;

    double eta = 0.05;     //!< driver rate
    double eps_rel = 0.01; //!< eps = eps_rel * max central amplitude ...
    double eps = -1.0;     //!< ... unless eps > 0 is given directly
    double kappa_fac = 0.5; //!< Omega smoothing kappa = kappa_fac * dx (level)
    double stiffness = 0.0; //!< d_t Omega += stiffness * lap Theta

    /// pure-gauge test: Theta(t=0) = a r^2/(r^2 + w^2) about a point on the
    /// symmetry axis (y = 0); physics must not change (use eta = kappa = 0)
    double test_theta0 = 0.0;
    double test_theta0_width = 5.0;
    double test_theta0_centre_x = 0.0;

    bool write_diagnostics = false; //!< write frame_diagnostics.dat (level 0)
};

#endif /* PHASEFRAMEPARAMS_HPP_ */
