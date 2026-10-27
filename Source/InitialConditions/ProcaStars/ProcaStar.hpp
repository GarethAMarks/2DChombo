/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PROCASTAR_HPP_
#define PROCASTAR_HPP_

#include "CCZ4CartoonProca.hpp"
#include "AsphericalPerturbation.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "ProcaStarBoost.hpp"
#include "ProcaStarSolution.hpp"
#include "MayDay.H"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "WeightFunction.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <array>
#include <cmath>
#include <complex>

//! Parameters of a Proca star, or of a head-on binary along the symmetry axis
struct ProcaStar_params_t
{
    ProcaStarSolution::params_t solver;          //!< 1D solver parameters (mu, A_central, eigen, perturbation, ...)
    std::array<double, CH_SPACEDIM> star_centre; //!< centre of the star (of the binary); must lie on the axis (y = 0)
    double phase;                                //!< phase of the field (of star 1), in units of pi
    bool write_1d_profiles;                      //!< write the 1D solution to ProcaStar_1d_profiles.dat (rank 0)
    double rapidity = 0.;                        //!< boost along x (> 0 moves towards -x)

    // binary (head-on along the axis), as the boson-star binaries of BosonStar.impl.hpp
    bool binary = false;
    ProcaStarSolution::params_t solver2;         //!< star 2 (defaults to star 1)
    double rapidity2 = 0.;                       //!< star 2 is boosted with -rapidity2 (> 0: towards +x)
    double phase2 = 0.;                          //!< phase of star 2, in units of pi
    double separation = 0.;                      //!< star 1 at x = q d/(1+q), star 2 at x = -d/(1+q)
    double mass_ratio = 1.;                      //!< q
    int id_choice = 2;                           //!< 0: plain superposition, 1: equal-mass fix, 2: unequal-mass fix
    double radius_width1 = 10., radius_width2 = 20.; //!< weight-function widths (id_choice 2)
    int conformal_factor_power = 4;              //!< used by id_choice 2

    //! aspherical (axisymmetric) perturbation of a single star, arXiv:2609.30021 Eq. (11), with phi the polar angle from
    //! the symmetry axis (AsphericalPerturbation.hpp); the spherical one is in solver (perturb, perturb_amp, ...)
    AsphericalPerturbation aspherical;
};

/// Initial data for a Proca star, or a head-on binary, in the 2D cartoon code
/**
 * Solves the spherical star with ProcaStarSolution (a port of the Spherical_BS_Evolver solver, optionally perturbed,
 * shared with ExoZvezda). For a star at rest, in isotropic coordinates centred on star_centre (on the symmetry axis),
 *   chi = psi^{-4}, h_ij = delta_ij, h_ww = 1, lapse = alpha, K = A_ij = Theta = shift = B = 0,
 *   X_phi = -f e^{i theta} / alpha,  X_i = i g n_i e^{i theta},  E^i = E n_i e^{i theta} / (alpha psi^4),  Z = 0,
 * with f, g, E the isotropic profiles of the 1D solver, n = (x, y)/r and theta = phase pi (|X_phi(0)| = A_central).
 * Boosts along the axis (rapidity) use the boson-star lab metric, lapse, shift and K_ij, and transform A_mu and F_mu nu
 * as a covector and a 2-form (ProcaStarBoost.hpp; derivation and checks in derivations/check_proca_boost.py). A boost
 * along the axis preserves axisymmetry (no azimuthal components). Binaries superpose X_phi, X_i and E^i and combine the
 * metrics as the boson-star binaries do (id_choice 0 plain, 1 equal-mass fix, 2 unequal-mass fix). Gamma^i must be
 * computed afterwards (GammaCartoonCalculator). Use with disable_simd() as for BosonStar.
 * Perturbations: the spherical one (m = 0) is applied by the 1D solver, with the constraints re-solved; the aspherical
 * one (aspherical.active, m >= 1; single stars only) multiplies the matter fields by exp(-2 delta cos(m phi)), phi the
 * polar angle from the symmetry axis about the star centre (arXiv:2609.30021, Eq. (11); AsphericalPerturbation.hpp).
 **/
class ProcaStar
{
  public:
    ProcaStar(ProcaStar_params_t a_params, double a_dx, int a_verbosity)
        : m_params(a_params), m_dx(a_dx), m_verbosity(a_verbosity)
    {
    }

    //! Computes the 1d solution(s) (every rank solves them; they are cheap)
    void compute_1d_solution()
    {
        if (std::abs(m_params.star_centre[1]) > 1e-12)
            MayDay::Error("ProcaStar: star_centre must lie on the symmetry axis (y = 0) in the cartoon code");
        if (!m_1d_sol.solve(m_params.solver, pout(), m_verbosity == 0))
            MayDay::Error("ProcaStar: the 1D Proca star solver failed (see pout); try a different "
                          "proca_frequency_guess / proca_alpha_central_guess or a larger proca_solver_radius");
        report(m_1d_sol, "ProcaStar");
        if (m_params.binary)
        {
            if (!m_1d_sol2.solve(m_params.solver2, pout(), m_verbosity == 0))
                MayDay::Error("ProcaStar: the 1D Proca star solver failed for star 2 (see pout)");
            report(m_1d_sol2, "ProcaStar 2");
        }
#ifdef CH_MPI
        int rank;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
        int rank = 0;
#endif
        if (m_params.write_1d_profiles && rank == 0)
            m_1d_sol.write_profiles("ProcaStar_1d_profiles.dat");
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        CCZ4CartoonProcaVars::VarsWithGauge<data_t> vars;
        VarsTools::assign(vars, 0.);
        Coordinates<data_t> coords(current_cell, m_dx, m_params.star_centre);
        const double x = coords.x; // data_t = double (disable_simd)
        const double y = coords.y;

        if (!m_params.binary)
        {
            const BoostedProcaStarVars b = boosted_proca_star(m_1d_sol, x, y, 0., m_params.rapidity, m_params.phase);
            add_field(vars, b, m_params.aspherical.factor(x, y, 0.)); // factor 1 unless perturbed aspherically
            set_metric(vars, b.gLL, b.KLL, b.lapse, b.shiftx, -1.);
        }
        else
        {
            const double q = m_params.mass_ratio, d = m_params.separation;
            const double xi1 = m_params.rapidity, xi2 = -m_params.rapidity2;
            const BoostedProcaStarVars s1 = boosted_proca_star(m_1d_sol, x - q * d / (q + 1.), y, 0., xi1, m_params.phase);
            const BoostedProcaStarVars s2 = boosted_proca_star(m_1d_sol2, x + d / (q + 1.), y, 0., xi2, m_params.phase2);
            add_field(vars, s1);
            add_field(vars, s2);

            double gLL[3][3], KLL[3][3], g1UU[3][3], g2UU[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                {
                    gLL[i][j] = s1.gLL[i][j] + s2.gLL[i][j] - ((i == j) ? 1. : 0.);
                    KLL[i][j] = 0.;
                }
            if (m_params.id_choice == 1)
            {
                // subtract star 1's metric at the position of star 2 (Helfer et al.)
                const BoostedProcaStarVars h = boosted_proca_star(m_1d_sol, -d, 0., 0., xi1, 0.);
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        gLL[i][j] = s1.gLL[i][j] + s2.gLL[i][j] - h.gLL[i][j];
            }
            inverse3(s1.gLL, g1UU);
            inverse3(s2.gLL, g2UU);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        for (int l = 0; l < 3; ++l)
                            KLL[i][j] += gLL[i][l] * (g1UU[l][k] * s1.KLL[k][j] + g2UU[l][k] * s2.KLL[k][j]);

            double chi_override = -1.;
            if (m_params.id_choice == 2)
            {
                // the unequal-mass fix of BosonStar.impl.hpp (id_choice 2): conformal factor corrected at both centres
                const BoostedProcaStarVars c1 = boosted_proca_star(m_1d_sol, -d, 0., 0., xi1, 0.);
                const BoostedProcaStarVars c2 = boosted_proca_star(m_1d_sol2, d, 0., 0., xi2, 0.);
                const BoostedProcaStarVars c11 = boosted_proca_star(m_1d_sol, 0., 0., 0., xi1, 0.);
                const BoostedProcaStarVars c22 = boosted_proca_star(m_1d_sol2, 0., 0., 0., xi2, 0.);
                const double n_power = m_params.conformal_factor_power / 12.0;
                auto det_sup = [](const double a[3][3], const double b[3][3]) {
                    return (a[0][0] + b[0][0] - 1.) * (a[1][1] + b[1][1] - 1.) * (a[2][2] + b[2][2] - 1.);
                };
                const double chi_plain = pow(gLL[0][0] * gLL[1][1] * gLL[2][2], n_power);
                const double delta_1 = pow(c11.gLL[0][0] * c11.gLL[1][1] * c11.gLL[2][2], n_power) -
                                       pow(det_sup(c11.gLL, c2.gLL), n_power);
                const double delta_2 = pow(c22.gLL[0][0] * c22.gLL[1][1] * c22.gLL[2][2], n_power) -
                                       pow(det_sup(c22.gLL, c1.gLL), n_power);
                WeightFunction weight;
                const double w1 = m_params.radius_width1, w2 = m_params.radius_width2;
                const double profile1 = weight.profile_chi((x - q * d / (q + 1)) * cosh(xi1), y, 0., w1);
                const double profile2 = weight.profile_chi((x + d / (q + 1)) * cosh(xi2), y, 0., w2);
                const double profile_11 = weight.profile_chi(0., 0., 0., w1);
                const double profile_12 = weight.profile_chi((d / (q + 1)) * (cosh(xi2) + q * cosh(xi1)), 0., 0., w1);
                const double profile_21 = weight.profile_chi((d / (q + 1)) * (-q * cosh(xi1) - cosh(xi2)), 0., 0., w2);
                const double profile_22 = weight.profile_chi(0., 0., 0., w2);
                const double det = profile_11 * profile_22 - profile_12 * profile_21;
                const double value1 = (-profile_21 * delta_2 + profile_22 * delta_1) / det;
                const double value2 = (profile_11 * delta_2 - profile_12 * delta_1) / det;
                const double chi_corrected = chi_plain + profile1 * value1 + profile2 * value2;
                chi_override = pow(chi_corrected, -4.0 / m_params.conformal_factor_power);
            }
            set_metric(vars, gLL, KLL, sqrt(s1.lapse * s1.lapse + s2.lapse * s2.lapse - 1.), s1.shiftx + s2.shiftx,
                       chi_override);
        }
        current_cell.store_vars(vars);
    }

    ProcaStarSolution m_1d_sol;
    ProcaStarSolution m_1d_sol2;

    //! adds the in-plane Proca variables of a (boosted) star at a point of the plane z = 0 (its X_z, E^z vanish there)
    //! (multiplied by F, the aspherical perturbation factor)
    template <class vars_t> static void add_field(vars_t &vars, const BoostedProcaStarVars &b, double F = 1.)
    {
        for (int i = 0; i < 2; ++i)
        {
            vars.X_Re[i] += F * b.X[i].real();
            vars.X_Im[i] += F * b.X[i].imag();
            vars.E_Re[i] += F * b.E[i].real();
            vars.E_Im[i] += F * b.E[i].imag();
        }
        vars.Xphi_Re += F * b.Xphi.real();
        vars.Xphi_Im += F * b.Xphi.imag();
    }

  protected:
    ProcaStar_params_t m_params;
    double m_dx;
    int m_verbosity;

    static void report(const ProcaStarSolution &s, const char *name)
    {
        pout() << name << ": A_central = " << s.A_central << ", f(0) = " << s.f0 << ", omega = " << s.omega << ", M = " << s.M
               << ", N = " << s.noether_charge << ", binding energy = " << s.binding_energy << ", r_99 = " << s.r_99
               << endl;
    }

    static void inverse3(const double g[3][3], double gi[3][3])
    {
        const double det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                           g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                           g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
        gi[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) / det;
        gi[0][1] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) / det;
        gi[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / det;
        gi[1][0] = (g[1][2] * g[2][0] - g[1][0] * g[2][2]) / det;
        gi[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) / det;
        gi[1][2] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) / det;
        gi[2][0] = (g[1][0] * g[2][1] - g[1][1] * g[2][0]) / det;
        gi[2][1] = (g[0][1] * g[2][0] - g[0][0] * g[2][1]) / det;
        gi[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) / det;
    }

    //! CCZ4 cartoon variables from gamma_ij, K_ij (3x3 with w = z), lapse and shift; chi_override > 0 replaces chi
    //! (the unequal-mass fix), and h_ij, A_ij are then built with it, as in BosonStar.impl.hpp
    template <class vars_t>
    static void set_metric(vars_t &vars, const double gLL[3][3], const double KLL[3][3], double lapse, double shiftx,
                           double chi_override)
    {
        double gUU[3][3];
        inverse3(gLL, gUU);
        const double chi_plain = pow(gLL[0][0] * gLL[1][1] * gLL[2][2], -1. / 3.);
        vars.chi = (chi_override > 0.) ? chi_override : chi_plain;
        const double chi_h = (chi_override > 0.) ? chi_override : chi_plain;
        double K = 0.;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                K += KLL[i][j] * gUU[i][j];
        vars.K = K;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
            {
                vars.h[i][j] = chi_h * gLL[i][j];
                vars.A[i][j] = chi_h * (KLL[i][j] - K * gLL[i][j] / 3.);
            }
        vars.hww = chi_h * gLL[2][2];
        vars.Aww = chi_h * (KLL[2][2] - K * gLL[2][2] / 3.);
        vars.lapse = lapse;
        vars.shift[0] = shiftx;
    }
};

#endif /* PROCASTAR_HPP_ */
