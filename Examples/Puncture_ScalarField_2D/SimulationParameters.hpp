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
//#include "ComplexPotential.hpp"

#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif

#ifdef USE_TWOPUNCTURES
#include "TP_Parameters.hpp"
#include "BoostedBH.hpp"
#endif


class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    #ifdef USE_TWOPUNCTURES
    double tp_offset_plus, tp_offset_minus;

    //param sets for TP data and each boosted BH
    TP::Parameters tp_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
    #endif


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
	pp.load("threshold_rho", regrid_threshold_rho, 1.);
        pp.load("threshold_chi", regrid_threshold_chi, 1.); 
        pp.load("regrid_wall_width", m_regrid_wall_width, 3.);
        pp.load("wall_min_regrid_level", m_wall_min_regrid_level, 3);
        pp.load("away_max_regrid_level", m_away_max_regrid_level, 0);

        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);

        // Puncture tracking tagging
        pp.load("do_puncture_track", do_puncture_track, false);
        pp.load("puncture1_level", puncture1_level, 6);
        pp.load("puncture2_level", puncture2_level, 6);
        puncture_max_levels[0] = puncture1_level;
        puncture_max_levels[1] = puncture2_level;
                                                      
        pp.load("tag_radius_A", tag_radius_A, 10.);
        pp.load("tag_radius_B", tag_radius_B, 10.);
        pp.load("tag_buffer", tag_buffer, 0.5);
        pp.load("puncture_min_separation", puncture_min_separation, 1e-3);
        //pp.load("use_TP_data", use_TP_data, false);

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


        //TwoPunctures parameters
        #ifdef USE_TWOPUNCTURES

	double v1 = -tanh(bosonstar_params.BS_rapidity);
	double v2 = tanh(bosonstar2_params.BS_rapidity);
	
	double gamma1 = 1 / sqrt ( 1 - v1 * v1);
	double gamma2 = 1 / sqrt ( 1 - v2 * v2);

    	tp_params.verbose = (verbosity > 0);	
	 
	double q = bosonstar_params.mass_ratio;
        double d = bosonstar_params.BS_separation;
	double b = bosonstar_params.BS_impact_parameter;
	
	bool calculate_target_masses;
        pp.load("TP_calculate_target_masses", calculate_target_masses, false);
        tp_params.give_bare_mass = !calculate_target_masses;
	
	// masses
        if (calculate_target_masses)
        {
            pp.load("TP_target_mass_plus", tp_params.target_M_plus, bosonstar_params.BlackHoleMass);
            pp.load("TP_target_mass_minus", tp_params.target_M_minus, bosonstar_params.BlackHoleMass2);
            pp.load("TP_adm_tol", tp_params.adm_tol, 1e-10);
            pout() << "The black holes have target ADM masses of "
                   << tp_params.target_M_plus << " and "
                   << tp_params.target_M_minus << "\n";
            bh2_params.mass = tp_params.target_M_minus;
            bh1_params.mass = tp_params.target_M_plus;
        }
        else
        {
            pp.load("TP_mass_plus", tp_params.par_m_plus);
            pp.load("TP_mass_minus", tp_params.par_m_minus);
            bh1_params.mass = tp_params.par_m_plus;
            bh2_params.mass = tp_params.par_m_minus;
            pout() << "The black holes have bare masses of "
                   << std::setprecision(16) << tp_params.par_m_plus << " and "
                   << tp_params.par_m_minus << "\n";
            // reset precision
            pout() << std::setprecision(6);
        }

	// BH spin and momenta
	std::array<double,2> momentum1{bosonstar_params.BlackHoleMass * gamma1 * v1, 0};
	std::array<double,2> momentum2{bosonstar_params.BlackHoleMass2 * gamma2 * v2, 0};

        //read in momenta and spins from input file, using zero for the z components as 2D
        std::array<double, CH_SPACEDIM> spin_minus, spin_plus;
        pp.load("TP_momentum_minus", bh2_params.momentum, momentum2);
        pp.load("TP_momentum_plus", bh1_params.momentum, momentum1);
        pp.load("TP_spin_plus", spin_plus);
        pp.load("TP_spin_minus", spin_minus);
        FOR(i)
        {
            tp_params.par_P_minus[i] = bh2_params.momentum[i];
            tp_params.par_P_plus[i] = bh1_params.momentum[i];
            tp_params.par_S_minus[i] = spin_minus[i];
            tp_params.par_S_plus[i] = spin_plus[i];
        }
        tp_params.par_P_minus[2] = 0.0;
        tp_params.par_P_plus[2] = 0.0;
        tp_params.par_S_minus[2] = 0.0;
        tp_params.par_S_plus[2] = 0.0;


        pout() << "The corresponding momenta are:";
        pout() << "\nP_plus = ";
        FOR(i) { pout() << tp_params.par_P_plus[i] << " "; }
        pout() << "\nP_minus = ";
        FOR(i) { pout() << tp_params.par_P_minus[i] << " "; }

        pout() << "\nThe corresponding spins are:";
        pout() << "\nS_plus = ";
        FOR(i) { pout() << tp_params.par_S_plus[i] << " "; }
        pout() << "\nS_minus = ";
        FOR(i) { pout() << tp_params.par_S_minus[i] << " "; }
        pout() << "\n";

	// interpolation type
        bool use_spectral_interpolation;
        pp.load("TP_use_spectral_interpolation", use_spectral_interpolation,
                false);
        tp_params.grid_setup_method =
            (use_spectral_interpolation) ? "evaluation" : "Taylor expansion";

        // initial_lapse (default to psi^n)
        pp.load("TP_initial_lapse", tp_params.initial_lapse,
                std::string("psi^n"));
        if (tp_params.initial_lapse != "twopunctures-antisymmetric" &&
            tp_params.initial_lapse != "twopunctures-averaged" &&
            tp_params.initial_lapse != "psi^n" &&
            tp_params.initial_lapse != "brownsville")
        {
            std::string message = "Parameter: TP_initial_lapse: ";
            message += tp_params.initial_lapse;
            message += " invalid";
            MayDay::Error(message.c_str());
        }
        if (tp_params.initial_lapse == "psi^n")
        {
            pp.load("TP_initial_lapse_psi_exponent",
                    tp_params.initial_lapse_psi_exponent, -2.0);
        }

        // Spectral grid parameters
        pp.load("TP_npoints_A", tp_params.npoints_A, 30);
        pp.load("TP_npoints_B", tp_params.npoints_B, 30);
        pp.load("TP_npoints_phi", tp_params.npoints_phi, 16);
        if (tp_params.npoints_phi % 4 != 0)
        {
            MayDay::Error("TP_npoints_phi must be a multiple of 4");
        }

        // Solver parameters and tolerances
        pp.load("TP_Newton_tol", tp_params.Newton_tol, 1e-10);
        pp.load("TP_Newton_maxit", tp_params.Newton_maxit, 5);
        pp.load("TP_epsilon", tp_params.TP_epsilon, 1e-6);
        pp.load("TP_Tiny", tp_params.TP_Tiny, 0.0);
        pp.load("TP_Extend_Radius", tp_params.TP_Extend_Radius, 0.0);

	//total distance between BH and BS
	double total_sep = sqrt(d*d + b*b);

        // BH positions
        pp.load("TP_offset_minus", tp_offset_minus,  -total_sep / (q + 1));
        pp.load("TP_offset_plus", tp_offset_plus, q * total_sep / (q + 1));
        bh1_params.center = center;
        bh2_params.center = center;
        bh1_params.center[0] += tp_offset_minus;
        bh2_params.center[0] += tp_offset_plus;

        //test: account for y offset due to symmetry enforcement
        //bh1_params.center[1] += L/2;
        //bh2_params.center[1] += L/2;


        double center_offset_x = 0.5 * (tp_offset_plus + tp_offset_minus);
        tp_params.center_offset[0] = center_offset_x;
        // par_b is half the distance between BH_minus and BH_plus
        tp_params.par_b = 0.5 * (tp_offset_plus - tp_offset_minus);
        pp.load("TP_swap_xz", tp_params.swap_xz, false);

        // Debug output
        pp.load("TP_do_residuum_debug_output",
                tp_params.do_residuum_debug_output, false);
        pp.load("TP_do_initial_debug_output", tp_params.do_initial_debug_output,
                false);

        // Irrelevant parameters set to default value
        tp_params.keep_u_around = false;
        tp_params.use_sources = false;
        tp_params.rescale_sources = true;
        tp_params.use_external_initial_guess = false;
        tp_params.multiply_old_lapse = false;
        tp_params.schedule_in_ADMBase_InitialData = true;
        tp_params.solve_momentum_constraint = false;
        tp_params.metric_type = "something else";
        tp_params.conformal_storage = "not conformal at all";
        tp_params.conformal_state = 0;
        tp_params.mp = 0;
        tp_params.mm = 0;
        tp_params.mp_adm = 0;
        tp_params.mm_adm = 0;
	#endif

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

    void check_params()
    {}

    // tagging
    bool activate_extraction;

    // For PhiAndK regridding
    double m_threshold_phi, m_threshold_K, m_threshold_rho, m_regrid_wall_width;
    int m_wall_min_regrid_level, m_away_max_regrid_level;

    // Layer regridding
    bool m_do_layer_tagging;
    double y_regrid_lim;

    bool do_puncture_track;
    int puncture1_level;
    int puncture2_level;
    std::array<int, 2> puncture_max_levels = {6, 6};
    double tag_radius_A;
    double tag_radius_B;
    double tag_buffer;
    double puncture_min_separation;

    //bool use_TP_data; //whether to use TwoPunctures initial data
    
    BosonStar_params_t bosonstar_params;
    BosonStar_params_t bosonstar2_params;
    Potential::params_t potential_params;

    double m_G_Newton;
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;
    bool identical;
    
    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_rho;

#ifdef USE_AHFINDER
    AHInitialGuessEllipsoid AH_1_initial_guess_ellipsoid;
    AHInitialGuessEllipsoid AH_2_initial_guess_ellipsoid;
    bool AH_set_origins_to_punctures;
#endif
};
#endif /* SIMULATIONPARAMETERS_HPP_ */
