/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "ProcaStar2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// Tagging
#include "ProcaAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "CCZ4CartoonProca.hpp"
#include "ConstraintsCartoonProca.hpp"
#include "GammaCartoonCalculator.hpp"
#include "MovingPunctureGauge.hpp"
#include "TraceARemovalCartoon.hpp"

// Proca star initial data
#include "ProcaStar.hpp"

// GW extraction and ADM quantities
#include "ADMQuantities.hpp"
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

void ProcaStar2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()), m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

void ProcaStar2DLevel::initialData()
{
    CH_TIME("ProcaStar2DLevel::initialData");
    if (m_verbosity)
        pout() << "ProcaStar2DLevel::initialData " << m_level << endl;

    ProcaStar proca_star(m_p.procastar_params, m_dx, m_verbosity);
    proca_star.compute_1d_solution();

    BoxLoops::loop(make_compute_pack(SetValue(0.0), proca_star), m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
}

// Things to do before a plot level - Weyl scalars, constraints and the Proca diagnostics
void ProcaStar2DLevel::prePlotLevel()
{
    fillAllGhosts();
    BoxLoops::loop(make_compute_pack(WeylOmScalar(m_p.extraction_params.center, m_dx),
                                     ConstraintsCartoonProca(m_dx, m_p.proca_params, m_p.m_G_Newton)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void ProcaStar2DLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs, const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()), a_soln, a_soln,
                   INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side with the Proca field
    BoxLoops::loop(CCZ4CartoonProca<MovingPunctureGauge, FourthOrderDerivatives>(
                       m_p.ccz4_params, m_dx, m_p.sigma, m_p.proca_params, m_p.m_G_Newton, m_p.formulation),
                   a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void ProcaStar2DLevel::specificUpdateODE(GRLevelData &a_soln, const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ProcaStar2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion, const FArrayBox &current_state,
                                               const double a_maximum)
{
    BoxLoops::loop(ProcaAndChiExtractionTaggingCriterion(m_dx, m_level, m_p.extraction_params,
                                                         m_p.regrid_threshold_Xphi, m_p.regrid_threshold_chi,
                                                         m_p.activate_extraction),
                   current_state, tagging_criterion);
}

void ProcaStar2DLevel::specificPostTimeStep()
{
    CH_TIME("ProcaStar2DLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.); // called during setup at t=0 from Main

    fillAllGhosts();

    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx), m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(ConstraintsCartoonProca(m_dx, m_p.proca_params, m_p.m_G_Newton), m_state_new,
                   m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    if (m_p.activate_extraction == 1 &&
        at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
                pout() << "ProcaStar2DLevel::specificPostTimeStep: Extracting gravitational waves." << endl;

            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time, first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    // integrands of the ADM surface integrals (plottable as M_adm, P_adm; to be integrated over an extraction sphere)
    if (at_level_timestep_multiple(0))
        BoxLoops::loop(ADMQuantities(m_p.extraction_params.center, m_dx, c_Madm, Interval(c_Padm, c_Padm)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    if (m_level == 0)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom2), 2, true);
        double L2_Gauss = amr_reductions.norm(c_Gauss, 2, true);
        SmallDataIO constraints_file("constraint_norms", m_dt, m_time, m_restart_time, SmallDataIO::APPEND,
                                     first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom", "L^2_Gauss"});
        constraints_file.write_time_data_line({L2_Ham, L2_Mom, L2_Gauss});

        // particle number (conserved up to outflow through the boundary)
        double particle_number = amr_reductions.sum(c_N);
        SmallDataIO N_file("ParticleNumber", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        N_file.remove_duplicate_time_data();
        if (first_step)
            N_file.write_header_line({"N"});
        N_file.write_time_data_line({particle_number});

        AMRReductions<VariableType::evolution> amr_reductions_evolution(m_gr_amr);
        double min_chi = amr_reductions_evolution.min(c_chi);
        SmallDataIO min_chi_file("min_chi", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (first_step)
            min_chi_file.write_header_line({"min_chi"});
        min_chi_file.write_time_data_line({min_chi});

        // maximum of |X_phi|, the central amplitude A_central = f(0)/alpha(0) of the 1D solver for a star at rest
        double mod_Xphi_max = amr_reductions.max(c_mod_Xphi);
        SmallDataIO mod_Xphi_file("mod_Xphi_max", m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step);
        mod_Xphi_file.remove_duplicate_time_data();
        if (first_step)
            mod_Xphi_file.write_header_line({"max mod Xphi"});
        mod_Xphi_file.write_time_data_line({mod_Xphi_max});
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif
}
