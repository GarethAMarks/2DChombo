/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "DiracStar2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// Tagging
#include "DiracAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "CCZ4CartoonDirac.hpp"
#include "ConstraintsCartoonDirac.hpp"
#include "GammaCartoonCalculator.hpp"
#include "MovingPunctureGauge.hpp"
#include "TraceARemovalCartoon.hpp"

// Dirac star initial data
#include "DiracStar.hpp"

// GW extraction and ADM quantities
#include "ADMQuantities.hpp"
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

void DiracStar2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()), m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

void DiracStar2DLevel::initialData()
{
    CH_TIME("DiracStar2DLevel::initialData");
    if (m_verbosity)
        pout() << "DiracStar2DLevel::initialData " << m_level << endl;

    DiracStar dirac_star(m_p.diracstar_params, m_dx, m_verbosity);
    dirac_star.compute_1d_solution();

    BoxLoops::loop(make_compute_pack(SetValue(0.0), dirac_star), m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
}

// Things to do before a plot level - Weyl scalars, constraints and the Dirac diagnostics
void DiracStar2DLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(make_compute_pack(WeylOmScalar(m_p.extraction_params.center, m_dx),
                                     ConstraintsCartoonDirac(m_dx, m_p.dirac_params, m_p.m_G_Newton)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void DiracStar2DLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs, const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()), a_soln, a_soln,
                   INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side with the Dirac field
    BoxLoops::loop(CCZ4CartoonDirac<MovingPunctureGauge, FourthOrderDerivatives>(
                       m_p.ccz4_params, m_dx, m_p.sigma, m_p.dirac_params, m_p.m_G_Newton, m_p.formulation),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void DiracStar2DLevel::specificUpdateODE(GRLevelData &a_soln, const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void DiracStar2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion, const FArrayBox &current_state,
                                               const double a_maximum)
{
    BoxLoops::loop(DiracAndChiExtractionTaggingCriterion(m_dx, m_level, m_p.extraction_params,
                                                         m_p.regrid_threshold_psi, m_p.regrid_threshold_chi,
                                                         m_p.activate_extraction),
                   current_state, tagging_criterion);
}

void DiracStar2DLevel::specificPostTimeStep()
{
    CH_TIME("DiracStar2DLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.); // called during setup at t=0 from Main

    fillAllGhosts();

    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx), m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(ConstraintsCartoonDirac(m_dx, m_p.dirac_params, m_p.m_G_Newton), m_state_new,
                   m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    if (m_p.activate_extraction == 1 &&
        at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
                pout() << "DiracStar2DLevel::specificPostTimeStep: Extracting gravitational waves." << endl;

            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time, first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    if (at_level_timestep_multiple(0))
        BoxLoops::loop(ADMQuantities(m_p.extraction_params.center, m_dx, 0), m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);

    if (m_level == 0)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom2), 2, true);
        SmallDataIO constraints_file("constraint_norms", m_dt, m_time, m_restart_time, SmallDataIO::APPEND,
                                     first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});

        double M_ADM = amr_reductions.sum(c_Madm);
        SmallDataIO M_ADM_file("M_ADM", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        M_ADM_file.remove_duplicate_time_data();
        if (first_step)
            M_ADM_file.write_header_line({"M_ADM"});
        M_ADM_file.write_time_data_line({M_ADM});

        // particle number of both spinors (conserved up to outflow through the boundary)
        double particle_number = amr_reductions.sum(c_N);
        SmallDataIO N_file("ParticleNumber", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        N_file.remove_duplicate_time_data();
        if (first_step)
            N_file.write_header_line({"N (both spinors)"});
        N_file.write_time_data_line({particle_number});

        AMRReductions<VariableType::evolution> amr_reductions_evolution(m_gr_amr);
        double min_chi = amr_reductions_evolution.min(c_chi);
        SmallDataIO min_chi_file("min_chi", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (first_step)
            min_chi_file.write_header_line({"min_chi"});
        min_chi_file.write_time_data_line({min_chi});

        // maximum of sqrt(4 pi psi^dag psi), the central amplitude A_central of the 1D solver for a star at rest
        double mod_psi_max = amr_reductions.max(c_mod_psi);
        SmallDataIO mod_psi_file("mod_psi_max", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        mod_psi_file.remove_duplicate_time_data();
        if (first_step)
            mod_psi_file.write_header_line({"max mod psi"});
        mod_psi_file.write_time_data_line({mod_psi_max});
    }
}
