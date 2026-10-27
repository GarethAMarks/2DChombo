/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "ProcaFieldCartoon.hpp"
#include "ProcaStar.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        pp.load("G_Newton", m_G_Newton, 1.0);

        // Proca field
        pp.load("proca_mass", proca_params.mass, 1.0);
        pp.load("proca_damping", proca_params.damping, 1.0); // kappa of the Gauss-constraint damping

        // Proca star initial data (parameter meanings as in Spherical_BS_Evolver's BSParams.par)
        ProcaStarSolution::params_t &s = procastar_params.solver;
        s.mu = proca_params.mass;
        s.G_Newton = m_G_Newton;
        pp.load("central_amplitude_proca", s.A_central);
        pp.load("proca_eigen", s.eigen, 0);
        pp.load("proca_lambda", s.lambda, 0.);
        pp.load("proca_solitonic", s.solitonic, false);
        pp.load("proca_sigma", s.sigma, 0.2);
        pp.load("proca_solver_radius", s.R, 150.);
        pp.load("proca_gridpoints", s.n_gridpoints, 30001);
        pp.load("proca_frequency_guess", s.frequency_guess, 1.5);
        pp.load("proca_alpha_central_guess", s.alpha_central, 0.9);
        pp.load("proca_r_match_fac", s.r_match_fac, 0.85);
        pp.load("proca_freq_epsilon", s.freq_epsilon, 1e-19);
        pp.load("proca_perturb", s.perturb, false);
        pp.load("proca_perturb_amp", s.perturb_amp, 0.);
        pp.load("proca_perturb_spread", s.perturb_spread, 3.);
        pp.load("proca_perturb_center", s.perturb_center, 0.);
        pp.load("star_centre", procastar_params.star_centre, {0.5 * L, 0.});
        pp.load("phase", procastar_params.phase, 0.0);
        pp.load("write_1d_profiles", procastar_params.write_1d_profiles, false);
        pp.load("PS_rapidity", procastar_params.rapidity, 0.0); // boost along the axis (> 0 moves towards -x)

        // head-on binary along the axis (conventions of the boson-star binaries in BosonStar.impl.hpp)
        pp.load("PS_binary", procastar_params.binary, false);
        procastar_params.solver2 = s;
        pp.load("central_amplitude_proca2", procastar_params.solver2.A_central, s.A_central);
        pp.load("PS_rapidity2", procastar_params.rapidity2, procastar_params.rapidity);
        pp.load("phase2", procastar_params.phase2, 0.0);
        pp.load("PS_separation", procastar_params.separation, 0.0);
        pp.load("mass_ratio", procastar_params.mass_ratio, 1.0);
        pp.load("id_choice", procastar_params.id_choice, 2);
        pp.load("radius_width1", procastar_params.radius_width1, 10.);
        pp.load("radius_width2", procastar_params.radius_width2, 20.);
        pp.load("conformal_factor_power", procastar_params.conformal_factor_power, 4);

        // tagging (as for ScalarField_2D: curvature of the matter field and of chi)
        pp.load("regrid_threshold_Xphi", regrid_threshold_Xphi, 1.e-4);
        pp.load("regrid_threshold_chi", regrid_threshold_chi, 1.);

        // Weyl extraction
        pp.load("activate_extraction", activate_extraction, false);
    }

    void check_params()
    {
        warn_parameter("proca_mass", proca_params.mass, proca_params.mass > 0., "must be positive");
        check_parameter("central_amplitude_proca", procastar_params.solver.A_central,
                        procastar_params.solver.A_central > 0., "must be positive");
        check_parameter("star_centre[1]", procastar_params.star_centre[1],
                        std::abs(procastar_params.star_centre[1]) < 1e-12,
                        "the star must be centred on the symmetry axis (y = 0)");
        check_parameter("proca_gridpoints", procastar_params.solver.n_gridpoints,
                        procastar_params.solver.n_gridpoints > 100, "must be > 100");
        // as in Spherical_BS_Evolver, only the free field is evolved
        check_parameter("proca_lambda", procastar_params.solver.lambda,
                        procastar_params.solver.lambda == 0. && !procastar_params.solver.solitonic,
                        "the evolution (ProcaFieldCartoon) is for the free Proca field: set proca_lambda = 0 and "
                        "proca_solitonic = 0");
        check_parameter("id_choice", procastar_params.id_choice,
                        procastar_params.id_choice >= 0 && procastar_params.id_choice <= 2, "must be 0, 1 or 2");
        warn_parameter("mass_ratio", procastar_params.mass_ratio,
                       !procastar_params.binary || procastar_params.id_choice != 1 || procastar_params.mass_ratio == 1.,
                       "the equal-mass fix (id_choice 1) assumes equal stars");
    }

    bool activate_extraction;

    double m_G_Newton;
    ProcaFieldCartoon::params_t proca_params;
    ProcaStar_params_t procastar_params;

    // Tagging thresholds
    Real regrid_threshold_Xphi, regrid_threshold_chi;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
