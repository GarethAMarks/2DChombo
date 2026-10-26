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
#include "DiracFieldCartoon.hpp"
#include "DiracStar.hpp"

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

        // Dirac field
        pp.load("dirac_mass", dirac_params.mass, 1.0);

        // Dirac star initial data (parameter meanings as in Spherical_BS_Evolver's BSParams.par)
        DiracStarSolution::params_t &s = diracstar_params.solver;
        s.mu = dirac_params.mass;
        s.G_Newton = m_G_Newton;
        pp.load("central_amplitude_dirac", s.A_central);
        pp.load("dirac_eigen", s.eigen, 0);
        pp.load("dirac_solver_radius", s.R, 100.);
        pp.load("dirac_gridpoints", s.n_gridpoints, 20001);
        pp.load("dirac_frequency_guess", s.frequency_guess, 1.8);
        pp.load("dirac_alpha_central_guess", s.alpha_central, 0.95);
        pp.load("dirac_r_match_fac", s.r_match_fac, 0.85);
        pp.load("dirac_freq_epsilon", s.freq_epsilon, 1e-19);
        pp.load("dirac_perturb", s.perturb, false);
        pp.load("dirac_perturb_amp", s.perturb_amp, 0.);
        pp.load("dirac_perturb_spread", s.perturb_spread, 3.);
        pp.load("dirac_perturb_center", s.perturb_center, 0.);
        pp.load("star_centre", diracstar_params.star_centre, {0.5 * L, 0.});
        pp.load("phase", diracstar_params.phase, 0.0);
        pp.load("write_1d_profiles", diracstar_params.write_1d_profiles, false);
        pp.load("DS_rapidity", diracstar_params.rapidity, 0.0); // boost along the axis (> 0 moves towards -x)

        // head-on binary along the axis (conventions of the boson-star binaries in BosonStar.impl.hpp)
        pp.load("DS_binary", diracstar_params.binary, false);
        diracstar_params.solver2 = s;
        pp.load("central_amplitude_dirac2", diracstar_params.solver2.A_central, s.A_central);
        pp.load("DS_rapidity2", diracstar_params.rapidity2, diracstar_params.rapidity);
        pp.load("phase2", diracstar_params.phase2, 0.0);
        pp.load("DS_separation", diracstar_params.separation, 0.0);
        pp.load("mass_ratio", diracstar_params.mass_ratio, 1.0);
        pp.load("id_choice", diracstar_params.id_choice, 2);
        pp.load("radius_width1", diracstar_params.radius_width1, 10.);
        pp.load("radius_width2", diracstar_params.radius_width2, 20.);
        pp.load("conformal_factor_power", diracstar_params.conformal_factor_power, 4);

        // tagging (as for ScalarField_2D: curvature of the matter field and of chi)
        pp.load("regrid_threshold_psi", regrid_threshold_psi, 1.e-4);
        pp.load("regrid_threshold_chi", regrid_threshold_chi, 1.);

        // Weyl extraction
        pp.load("activate_extraction", activate_extraction, false);
    }

    void check_params()
    {
        warn_parameter("dirac_mass", dirac_params.mass, dirac_params.mass > 0., "must be positive");
        check_parameter("central_amplitude_dirac", diracstar_params.solver.A_central,
                        diracstar_params.solver.A_central > 0., "must be positive");
        check_parameter("star_centre[1]", diracstar_params.star_centre[1],
                        std::abs(diracstar_params.star_centre[1]) < 1e-12,
                        "the star must be centred on the symmetry axis (y = 0)");
        check_parameter("dirac_gridpoints", diracstar_params.solver.n_gridpoints,
                        diracstar_params.solver.n_gridpoints > 100, "must be > 100");
        check_parameter("id_choice", diracstar_params.id_choice,
                        diracstar_params.id_choice >= 0 && diracstar_params.id_choice <= 2, "must be 0, 1 or 2");
        warn_parameter("mass_ratio", diracstar_params.mass_ratio,
                       !diracstar_params.binary || diracstar_params.id_choice != 1 || diracstar_params.mass_ratio == 1.,
                       "the equal-mass fix (id_choice 1) assumes equal stars");
    }

    bool activate_extraction;

    double m_G_Newton;
    DiracFieldCartoon::params_t dirac_params;
    DiracStar_params_t diracstar_params;

    // Tagging thresholds
    Real regrid_threshold_psi, regrid_threshold_chi;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
