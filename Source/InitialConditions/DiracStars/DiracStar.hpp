/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIRACSTAR_HPP_
#define DIRACSTAR_HPP_

#include "CCZ4CartoonDirac.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DiracStarBoost.hpp"
#include "DiracStarSolution.hpp"
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

//! Parameters of a Dirac star, or of a head-on binary along the symmetry axis
struct DiracStar_params_t
{
    DiracStarSolution::params_t solver;          //!< 1D solver parameters (mu, A_central, eigen, perturbation, ...)
    std::array<double, CH_SPACEDIM> star_centre; //!< centre of the star (of the binary); must lie on the axis (y = 0)
    double phase;                                //!< phase of the spinor (of star 1), in units of pi
    bool write_1d_profiles;                      //!< write the 1D solution to DiracStar_1d_profiles.dat (rank 0)
    double rapidity = 0.;                        //!< boost along x (> 0 moves towards -x)

    // binary (head-on along the axis), as the boson-star binaries of BosonStar.impl.hpp
    bool binary = false;
    DiracStarSolution::params_t solver2;         //!< star 2 (defaults to star 1)
    double rapidity2 = 0.;                       //!< star 2 is boosted with -rapidity2 (> 0: towards +x)
    double phase2 = 0.;                          //!< phase of star 2, in units of pi
    double separation = 0.;                      //!< star 1 at x = q d/(1+q), star 2 at x = -d/(1+q)
    double mass_ratio = 1.;                      //!< q
    int id_choice = 2;                           //!< 0: plain superposition, 1: equal-mass fix, 2: unequal-mass fix
    double radius_width1 = 10., radius_width2 = 20.; //!< weight-function widths (id_choice 2)
    int conformal_factor_power = 4;              //!< used by id_choice 2
};

/// Initial data for a Dirac star, or a head-on binary, in the 2D cartoon code
/**
 * Solves the spherical star with DiracStarSolution (a port of the Spherical_BS_Evolver solver, optionally perturbed)
 * and sets the J_x = +1/2 spinor (the J_x = -1/2 one is its mirror image, see DiracFieldCartoon.hpp). For a star at
 * rest, in isotropic coordinates centred on star_centre (on the symmetry axis),
 *   chi = psi^{-4}, h_ij = delta_ij, h_ww = 1, lapse = alpha, K = A_ij = Theta = shift = B = 0,
 *   psi0 = g e^{i theta} / sqrt(4 pi),  psi1 = 0,
 *   psi2 = -i f (x/r) e^{i theta} / sqrt(4 pi),  psi3 = -i f (y/r) e^{i theta} / sqrt(4 pi),   theta = phase pi.
 * Boosts along the axis (rapidity) use the boson-star lab metric, lapse, shift and K_ij (DiracStarBoost.hpp) and
 * the local Lorentz transformation of the spinor, psi_lab = U(eta) psi_rest(t~, x~), U = cosh(eta/2) - sinh(eta/2)
 * alpha_x with alpha_x block off-diagonal sigma_3 in this basis (derivation and checks in ExoZvezda's
 * Examples/SingleDiracStar/derivations/check_dirac_boost.py). A boost along the axis preserves axisymmetry and the
 * mirror relation between the two spinors. Binaries superpose the spinors and combine the metrics as the boson-star
 * binaries do (id_choice 0 plain, 1 equal-mass fix, 2 unequal-mass fix). Gamma^i must be computed afterwards
 * (GammaCartoonCalculator). Use with disable_simd() as for BosonStar.
 **/
class DiracStar
{
  public:
    DiracStar(DiracStar_params_t a_params, double a_dx, int a_verbosity)
        : m_params(a_params), m_dx(a_dx), m_verbosity(a_verbosity)
    {
    }

    //! Computes the 1d solution(s) (every rank solves them; they are cheap)
    void compute_1d_solution()
    {
        if (std::abs(m_params.star_centre[1]) > 1e-12)
            MayDay::Error("DiracStar: star_centre must lie on the symmetry axis (y = 0) in the cartoon code");
        if (!m_1d_sol.solve(m_params.solver, pout(), m_verbosity == 0))
            MayDay::Error("DiracStar: the 1D Dirac star solver failed (see pout); try a different "
                          "dirac_frequency_guess / dirac_alpha_central_guess or a larger dirac_solver_radius");
        report(m_1d_sol, "DiracStar");
        if (m_params.binary)
        {
            if (!m_1d_sol2.solve(m_params.solver2, pout(), m_verbosity == 0))
                MayDay::Error("DiracStar: the 1D Dirac star solver failed for star 2 (see pout)");
            report(m_1d_sol2, "DiracStar 2");
        }
#ifdef CH_MPI
        int rank;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
        int rank = 0;
#endif
        if (m_params.write_1d_profiles && rank == 0)
            m_1d_sol.write_profiles("DiracStar_1d_profiles.dat");
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        CCZ4CartoonDiracVars::VarsWithGauge<data_t> vars;
        VarsTools::assign(vars, 0.);
        Coordinates<data_t> coords(current_cell, m_dx, m_params.star_centre);
        const double x = coords.x; // data_t = double (disable_simd)
        const double y = coords.y;

        if (!m_params.binary)
        {
            const BoostedDiracStarVars b = boosted_dirac_star(m_1d_sol, x, y, 0., m_params.rapidity, m_params.phase);
            double psi_Re[4], psi_Im[4];
            spinor(b, psi_Re, psi_Im);
            for (int c = 0; c < 4; ++c)
            {
                vars.psi_Re[c] = psi_Re[c];
                vars.psi_Im[c] = psi_Im[c];
            }
            set_metric(vars, b.gLL, b.KLL, b.lapse, b.shiftx, -1.);
        }
        else
        {
            const double q = m_params.mass_ratio, d = m_params.separation;
            const double xi1 = m_params.rapidity, xi2 = -m_params.rapidity2;
            const BoostedDiracStarVars s1 = boosted_dirac_star(m_1d_sol, x - q * d / (q + 1.), y, 0., xi1, m_params.phase);
            const BoostedDiracStarVars s2 = boosted_dirac_star(m_1d_sol2, x + d / (q + 1.), y, 0., xi2, m_params.phase2);
            double p1R[4], p1I[4], p2R[4], p2I[4];
            spinor(s1, p1R, p1I);
            spinor(s2, p2R, p2I);
            for (int c = 0; c < 4; ++c)
            {
                vars.psi_Re[c] = p1R[c] + p2R[c];
                vars.psi_Im[c] = p1I[c] + p2I[c];
            }

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
                const BoostedDiracStarVars h = boosted_dirac_star(m_1d_sol, -d, 0., 0., xi1, 0.);
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
                const BoostedDiracStarVars c1 = boosted_dirac_star(m_1d_sol, -d, 0., 0., xi1, 0.);
                const BoostedDiracStarVars c2 = boosted_dirac_star(m_1d_sol2, d, 0., 0., xi2, 0.);
                const BoostedDiracStarVars c11 = boosted_dirac_star(m_1d_sol, 0., 0., 0., xi1, 0.);
                const BoostedDiracStarVars c22 = boosted_dirac_star(m_1d_sol2, 0., 0., 0., xi2, 0.);
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

    DiracStarSolution m_1d_sol;
    DiracStarSolution m_1d_sol2;

    //! the J_x = +1/2 spinor of a (boosted) star in the cartoon basis (spin along x: sigma'_x = sigma_3,
    //! sigma'_y = sigma_1), at a point of the plane z = 0
    static void spinor(const BoostedDiracStarVars &b, double *psi_Re, double *psi_Im)
    {
        using cd = std::complex<double>;
        const cd I(0., 1.);
        const cd pref = b.phase / sqrt(4. * M_PI);
        // rest frame: (g, 0, -i f n_x, -i f n_y)
        const cd u0 = b.g, u1 = 0., l0 = -I * b.f * b.n[0], l1 = -I * b.f * b.n[1];
        // U = cosh(eta/2) - sinh(eta/2) alpha_x, alpha_x = [[0, sigma_3], [sigma_3, 0]]
        const double ch = b.cosh_half_eta, sh = b.sinh_half_eta;
        const cd p[4] = {pref * (ch * u0 - sh * l0), pref * (ch * u1 + sh * l1), pref * (ch * l0 - sh * u0),
                         pref * (ch * l1 + sh * u1)};
        for (int c = 0; c < 4; ++c)
        {
            psi_Re[c] = p[c].real();
            psi_Im[c] = p[c].imag();
        }
    }

  protected:
    DiracStar_params_t m_params;
    double m_dx;
    int m_verbosity;

    static void report(const DiracStarSolution &s, const char *name)
    {
        pout() << name << ": A_central = " << s.A_central << ", omega = " << s.omega << ", M = " << s.M
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

#endif /* DIRACSTAR_HPP_ */
