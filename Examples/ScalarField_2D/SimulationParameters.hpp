/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#undef USE_AHFINDER //just ignore ah finder for now

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "ComplexPotential.hpp"

// Problem specific(from robin) includes:
#include "BosonStarParams.hpp"

// U(1) phase frame for the scalar field
#include "BoundaryConditions.hpp"
#include "PhaseFrameParams.hpp"
#include <algorithm>
#include <cmath>
//#include "ComplexPotential.hpp"

#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif

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
        pp.load("G_Newton", m_G_Newton);

        // Boson Star initial data params
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("eigen", bosonstar_params.eigen, 0);
        pp.load("gridpoints", bosonstar_params.gridpoints, 400000);
         pp.load("star_centre", bosonstar_params.star_centre,
                 {0.5 * L, 0.});

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_soliton", potential_params.sigma_soliton, 0.02);
        pp.load("BS_binary", bosonstar_params.BS_binary, false);
        pp.load("BS_BH_binary", bosonstar_params.BS_BH_binary, false);
        pp.load("BH_binary", bosonstar_params.BH_binary, false);
        pp.load("read_TS_1", bosonstar_params.read_TS_1, false);
        pp.load("read_TS_2", bosonstar_params.read_TS_2, false);
        pp.load("antiboson", bosonstar_params.antiboson, false);
        pp.load("BlackHoleMass", bosonstar_params.BlackHoleMass, 0.);
        pp.load("BlackHoleMass2", bosonstar_params.BlackHoleMass2, 0.);
        pp.load("BS_rapidity", bosonstar_params.BS_rapidity, 0.0);
         pp.load("BS_separation", bosonstar_params.BS_separation, 0.0);
        pp.load("BS_impact_parameter", bosonstar_params.BS_impact_parameter,
                0.0);
        pp.load("id_choice", bosonstar_params.id_choice, 2);
        pp.load("mass_ratio", bosonstar_params.mass_ratio, 1.0);
        pp.load("radius_width1", bosonstar_params.radius_width1, 10.);
        pp.load("radius_width2", bosonstar_params.radius_width2, 20.);
        pp.load("conformal_factor_power", bosonstar_params.conformal_factor_power, -4);
        pp.load("G_Newton", bosonstar_params.Newtons_constant, 1.0);
        pp.load("print_asymptotics", bosonstar_params.print_asymptotics, false);

        // Initialize values for bosonstar2_params to same as bosonstar_params
        // and then assign that ones that should differ below
        bosonstar2_params = bosonstar_params;

        // Are the two stars' profiles identical
        pp.load("identical", identical, false);

        // Boson Star 2 parameters
        if (!identical)
        {
            pp.load("central_amplitude_CSF2",
                    bosonstar2_params.central_amplitude_CSF);
            pp.load("BS_rapidity2",
                    bosonstar2_params.BS_rapidity);
        }

        // BubbleTaggingCriterion regridding
        pp.load("threshold_phi", regrid_threshold_phi, 1.);
        pp.load("threshold_chi", regrid_threshold_chi, 1.); 
        pp.load("regrid_wall_width", m_regrid_wall_width, 3.);
        pp.load("wall_min_regrid_level", m_wall_min_regrid_level, 3);
        pp.load("away_max_regrid_level", m_away_max_regrid_level, 0);

        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.0, 0.0});

        read_phase_frame_params(pp);

        // Charge projection of the KO dissipation of the scalar field
        pp.load("ko_charge_projection", ko_charge_projection, false);
        pp.load("ko_charge_projection_chi", ko_charge_projection_chi, true);

        // 4th-order (axis-corrected) Noether-charge weights
        pp.load("axis_corrected_charge", axis_corrected_charge, false);


#ifdef USE_AHFINDER
        pp.load("AH_set_origins_to_punctures", AH_set_origins_to_punctures,
                false);

        double guess1, guess2;
        pp.load("AH_1_initial_guess", guess1, 0.5 * bh1_params.mass);
        pp.load("AH_2_initial_guess", guess2, 0.5 * bh2_params.mass);

        double r_x_1 = guess1, r_y_1 = guess1;
        double r_x_2 = guess2, r_y_2 = guess2;

        double vel1 = bh1_params.momentum[0], vel2 = bh2_params.momentum[0];
        double contraction1 = sqrt(1. - vel1 * vel1),
               contraction2 = sqrt(1. - vel2 * vel2);

        double ah1_ellipsoid_contraction, ah2_ellipsoid_contraction;
        pp.load("AH_1_ellipsoid_contraction", ah1_ellipsoid_contraction,
                contraction1);
        pp.load("AH_2_ellipsoid_contraction", ah2_ellipsoid_contraction,
                contraction2);

        r_x_1 *= ah1_ellipsoid_contraction;
        r_x_2 *= ah2_ellipsoid_contraction;

        AH_1_initial_guess_ellipsoid.set_params(r_x_1, r_y_1);
        AH_2_initial_guess_ellipsoid.set_params(r_x_2, r_y_2);
#endif
    }

    /// U(1) phase frame parameters (see params.txt for descriptions)
    void read_phase_frame_params(GRParmParse &pp)
    {
        PhaseFrame_params_t &fp = phase_frame_params;
        pp.load("use_phase_frame", fp.use_phase_frame, false);
        pp.load("frame_links", fp.links, true);
        pp.load("frame_local_driver", fp.local_driver, true);
        pp.load("frame_omega_from_ID", fp.omega_from_ID, true);
        pp.load("frame_omega0", fp.omega0, 0.0);
        pp.load("frame_eta", fp.eta, 0.05);
        pp.load("frame_eps_rel", fp.eps_rel, 0.01);
        pp.load("frame_eps", fp.eps, -1.0);
        pp.load("frame_kappa_fac", fp.kappa_fac, 0.5);
        pp.load("frame_stiffness", fp.stiffness, 0.0);
        pp.load("frame_test_theta0", fp.test_theta0, 0.0);
        pp.load("frame_test_theta0_width", fp.test_theta0_width, 5.0);
        pp.load("frame_test_theta0_centre_x", fp.test_theta0_centre_x,
                bosonstar_params.star_centre[0]);
        pp.load("write_frame_diagnostics", fp.write_diagnostics, false);

        // eps: given directly, or relative to the largest central amplitude
        if (fp.eps <= 0.)
        {
            double amp_max = std::abs(bosonstar_params.central_amplitude_CSF);
            if (bosonstar_params.BS_binary)
                amp_max = std::max(
                    amp_max, std::abs(bosonstar2_params.central_amplitude_CSF));
            fp.eps = fp.eps_rel * amp_max;
        }
    }

    void check_params()
    {
        warn_parameter("ko_charge_projection", ko_charge_projection,
                       !ko_charge_projection || sigma != 0.,
                       "sigma = 0: no KO dissipation to project");

        const PhaseFrame_params_t &fp = phase_frame_params;
        if (fp.use_phase_frame)
        {
            check_parameter("frame_eps", fp.eps, fp.eps > 0.,
                            "must be > 0 (set frame_eps, or frame_eps_rel > 0 "
                            "with a nonzero central amplitude)");
            check_parameter("frame_eta", fp.eta, fp.eta >= 0., "must be >= 0");
            // explicit diffusion of Omega: kappa dt / dx^2 =
            // frame_kappa_fac * dt_multiplier must stay below ~0.26 (RK4,
            // 4th-order 2D Laplacian), less near the axis (d_y / y term)
            check_parameter("frame_kappa_fac", fp.kappa_fac,
                            fp.kappa_fac >= 0. &&
                                fp.kappa_fac * dt_multiplier <= 0.2,
                            "need 0 <= frame_kappa_fac * dt_multiplier <= 0.2 "
                            "for stable explicit smoothing of Omega");
            // relaxation term: eta * dt <= ~2.8 on the coarsest level
            warn_parameter("frame_eta", fp.eta,
                           fp.eta * coarsest_dt <= 2.,
                           "eta * coarsest dt > 2: driver may be unstable on "
                           "the coarsest level");
            check_parameter("frame_stiffness", fp.stiffness,
                            fp.stiffness >= 0., "must be >= 0");
            warn_parameter("frame_stiffness", fp.stiffness,
                           fp.stiffness <= 1.,
                           "Theta becomes a wave with speed "
                           "sqrt(frame_stiffness) > 1");
            warn_parameter("frame_test_theta0", fp.test_theta0,
                           fp.test_theta0 == 0.,
                           "pure-gauge test mode: results must agree with "
                           "use_phase_frame = 0 up to truncation error");
            warn_parameter("frame_local_driver", fp.local_driver,
                           fp.local_driver || !fp.omega_from_ID ||
                               fp.test_theta0 != 0.,
                           "fixed frame with Omega frozen at the local "
                           "initial-data frequency profile");

            // frame_theta and frame_omega should be extrapolated at outer
            // boundaries (MIXED_BC with extrapolating_vars), not Sommerfeld
            bool sommerfeld_bc = false;
            for (int idir = 0; idir < CH_SPACEDIM; ++idir)
            {
                if (!boundary_params.is_periodic[idir] &&
                    (boundary_params.hi_boundary[idir] ==
                         BoundaryConditions::SOMMERFELD_BC ||
                     boundary_params.lo_boundary[idir] ==
                         BoundaryConditions::SOMMERFELD_BC))
                    sommerfeld_bc = true;
            }
            warn_parameter("hi_boundary/lo_boundary", sommerfeld_bc,
                           !sommerfeld_bc,
                           "with use_phase_frame, use MIXED_BC (4) with "
                           "extrapolating_vars = frame_theta frame_omega "
                           "instead of Sommerfeld (1)");
        }
    }

    PhaseFrame_params_t phase_frame_params;
    bool ko_charge_projection;     //!< remove the charge-changing part of KO
    bool ko_charge_projection_chi; //!< include KO of chi (charge weight)
    bool axis_corrected_charge; //!< Euler-Maclaurin axis correction of N

    // tagging
    bool activate_extraction;

    // For PhiAndK regridding
    double m_threshold_phi, m_threshold_K, m_threshold_rho, m_regrid_wall_width;
    int m_wall_min_regrid_level, m_away_max_regrid_level;

    // Layer regridding
    bool m_do_layer_tagging;
    double y_regrid_lim;
    
    BosonStar_params_t bosonstar_params;
    BosonStar_params_t bosonstar2_params;
    Potential::params_t potential_params;

    double m_G_Newton;
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;
    bool identical;
    
    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi;

#ifdef USE_AHFINDER
    AHInitialGuessEllipsoid AH_1_initial_guess_ellipsoid;
    AHInitialGuessEllipsoid AH_2_initial_guess_ellipsoid;
    bool AH_set_origins_to_punctures;
#endif
};
#endif /* SIMULATIONPARAMETERS_HPP_ */
