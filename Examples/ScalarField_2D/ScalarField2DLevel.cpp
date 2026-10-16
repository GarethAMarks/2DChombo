/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include <fstream>
#include <mutex>
#include <memory>

#include "ScalarField2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "CCZ4Cartoon.hpp"
#include "MovingPunctureGauge.hpp"
#include "ConstraintsCartoon.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"
#include "NoetherCharge.hpp"

// For GW extraction
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

//Boson Star
#include "BosonStar.hpp"
#include "ComplexScalarField.hpp"
#include "ComplexPotential.hpp"
#include "ComputeWeightFunction.hpp"

#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"
#include "GammaCartoonCalculator.hpp"

// U(1) phase frame
#include "CCZ4CartoonPhaseFrame.hpp"
#include "PhaseFrameTools.hpp"

// Charge projection of the KO dissipation
#include "KOChargeProjection.hpp"
#include "BoxIterator.H"
#include "SPMD.H"
#include <vector>

void ScalarField2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void ScalarField2DLevel::initialData()
{
    CH_TIME("ScalarField2DLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarField2DLevel::initialData " << m_level << endl;

    // When changing class here, don't forget to change potential if necessary
    double spacing = .01;
    BosonStar boson_star(m_p.bosonstar_params, m_p.bosonstar2_params, m_p.potential_params,
                         m_p.m_G_Newton, m_dx, m_verbosity);


    if (m_verbosity)
        pout() << "BosonStarLevel::initialData - Compute_1d_solution " << m_level << endl;
    // the max radius the code might need to calculate out to is L*sqrt(3)
    boson_star.compute_1d_solution(4. * m_p.L);

    BoxLoops::loop(make_compute_pack(SetValue(0.0), boson_star), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());

    // Phase frame: frame_theta = 0 (psi = phi), frame_omega from the initial
    // data (in place: writes only the frame variables; the ghosts of
    // frame_omega are filled by fillAllGhosts below). For the pure-gauge test,
    // then set Theta_0(x) and rotate the field into that frame (pointwise).
    if (m_p.phase_frame_params.use_phase_frame)
    {
        BoxLoops::loop(PhaseFrameInitialOmega(m_p.phase_frame_params, m_dx),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
        if (m_p.phase_frame_params.test_theta0 != 0.)
            BoxLoops::loop(
                PhaseFrameSetTestTheta(m_p.phase_frame_params, m_dx),
                m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
    }

    BoxLoops::loop(ComputeWeightFunction(m_p.bosonstar_params, m_p.bosonstar2_params, m_dx), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
    
}

// Things to do before a plot level - need to calculate the Weyl scalars
void ScalarField2DLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    // diagnostics need the physical scalar field (identical to m_state_new
    // when use_phase_frame = 0)
    const GRLevelData &phys_state = physical_state();
    BoxLoops::loop(make_compute_pack(WeylOmScalar(m_p.extraction_params.center, m_dx), 
                                    Constraints<Potential>(m_dx, potential, m_p.m_G_Newton), 
                                   /* ADMQuantities(m_p.extraction_params.center, m_dx, 0),*/
                                    NoetherCharge<FourthOrderDerivatives>(m_dx, m_p.axis_corrected_charge),
                                    PhaseFrameCopyPhysical()), phys_state, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void ScalarField2DLevel::specificEvalRHS(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    Potential potential(m_p.potential_params);
    const PhaseFrame_params_t &frame_params = m_p.phase_frame_params;
    // the state the RHS (and its KO dissipation) was evaluated on
    const GRLevelData *rhs_input = &a_soln;
    GRLevelData phys_soln;
    if (!frame_params.use_phase_frame)
    {
        // original evolution; the frame variables stay at zero
        a_rhs.setVal(0., Interval(c_frame_theta, c_frame_omega));
        BoxLoops::loop(
            CCZ4Cartoon<MovingPunctureGauge, FourthOrderDerivatives, Potential>(
                m_p.ccz4_params, m_dx, m_p.sigma, potential, m_p.m_G_Newton,
                m_p.formulation),
            a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else
    {
        CCZ4CartoonPhaseFrame<MovingPunctureGauge, FourthOrderDerivatives,
                              Potential>
            frame_rhs(m_p.ccz4_params, m_dx, m_p.sigma, potential,
                      m_p.m_G_Newton, frame_params, m_p.formulation);
        if (frame_params.links)
        {
            // Link variables: give the RHS the physical field
            // exp(i Theta) psi in every cell (incl. ghosts, which were filled
            // with psi and Theta by the core), so all stencils act on the
            // physical field; CCZ4CartoonPhaseFrame rotates the result back
            // into the local frame.
            frame_to_physical(a_soln, phys_soln);
            BoxLoops::loop(frame_rhs, phys_soln, a_rhs, EXCLUDE_GHOST_CELLS);
            rhs_input = &phys_soln;
        }
        else
        {
            // continuum connection: works on psi directly
            BoxLoops::loop(frame_rhs, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
        }
    }

    // Remove the charge-changing part of the KO dissipation of the matter
    if (m_p.ko_charge_projection && m_p.sigma != 0.)
        ko_charge_projection(*rhs_input, a_soln, a_rhs);
}

void ScalarField2DLevel::ko_charge_projection(const GRLevelData &a_rhs_input,
                                              const GRLevelData &a_soln,
                                              GRLevelData &a_rhs)
{
    CH_TIME("ScalarField2DLevel::ko_charge_projection");

    // per-cell charge rate of the dissipation and normalisation (see
    // KOChargeProjection.hpp), from the field the dissipation acted on
    GRLevelData dens;
    dens.define(a_soln.disjointBoxLayout(), 2, IntVect::Zero);
    BoxLoops::loop(KOChargeRateDensity(m_dx, m_p.sigma,
                                       m_p.ko_charge_projection_chi,
                                       m_p.axis_corrected_charge),
                   a_rhs_input, dens, EXCLUDE_GHOST_CELLS);

    // only cells not covered by the next finer level count in the composite
    // charge
    std::vector<Box> covered;
    if (m_finer_level_ptr != nullptr)
    {
        const Vector<Box> &fine_boxes = m_finer_level_ptr->boxes();
        for (int i = 0; i < fine_boxes.size(); ++i)
        {
            Box b = fine_boxes[i];
            b.coarsen(m_ref_ratio);
            covered.push_back(b);
        }
    }

    double sums[2] = {0.0, 0.0};
    for (DataIterator dit = dens.dataIterator(); dit.ok(); ++dit)
    {
        const Box &valid = dens.disjointBoxLayout()[dit()];
        FArrayBox mask(valid, 1);
        mask.setVal(1.0);
        for (const Box &c : covered)
        {
            const Box overlap = valid & c;
            if (!overlap.isEmpty())
                mask.setVal(0.0, overlap, 0, 1);
        }
        const FArrayBox &d = dens[dit()];
        for (BoxIterator bit(valid); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            sums[0] += mask(iv, 0) * d(iv, 0);
            sums[1] += mask(iv, 0) * d(iv, 1);
        }
    }
#ifdef CH_MPI
    MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, Chombo_MPI::comm);
#endif
    if (!(sums[1] > 0.0))
        return; // no field on the uncovered part of this level

    // rhs -= lambda g,  g = (-Pi_Im, Pi, phi_Im, -phi) of the evolved state
    const double lambda = sums[0] / sums[1];
    for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
        const Box &valid = a_rhs.disjointBoxLayout()[dit()];
        FArrayBox &r = a_rhs[dit()];
        const FArrayBox &u = a_soln[dit()];
        for (BoxIterator bit(valid); bit.ok(); ++bit)
        {
            const IntVect iv = bit();
            r(iv, c_phi) += lambda * u(iv, c_Pi_Im);
            r(iv, c_phi_Im) -= lambda * u(iv, c_Pi);
            r(iv, c_Pi) -= lambda * u(iv, c_phi_Im);
            r(iv, c_Pi_Im) += lambda * u(iv, c_phi);
        }
    }

    if (m_verbosity)
        pout() << "ScalarField2DLevel::ko_charge_projection: level " << m_level
               << ", KO charge rate (uncovered, unnormalised) " << sums[0]
               << ", lambda " << lambda << endl;
}

void ScalarField2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarField2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state,
                                                 const double a_maximum)
{
    ComplexPhiAndChiExtractionTaggingCriterion criterion(m_dx, m_level,
                   m_p.mass_extraction_params, m_p.regrid_threshold_phi,
                   m_p.regrid_threshold_chi);
    if (!m_p.phase_frame_params.use_phase_frame)
    {
        BoxLoops::loop(criterion, current_state, tagging_criterion);
    }
    else
    {
        // the criterion uses derivatives of phi and Pi, which are not U(1)
        // invariant: tag on the physical field so refinement is unchanged
        FArrayBox phys_state(current_state.box(), current_state.nComp());
        phys_state.copy(current_state);
        BoxLoops::loop(PhaseFrameRotate(+1), phys_state, phys_state);
        BoxLoops::loop(criterion, phys_state, tagging_criterion);
    }
}

void ScalarField2DLevel::specificPostTimeStep()
{
    CH_TIME("ScalarField2DLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main

    fillAllGhosts();

    // diagnostics need the physical scalar field (identical to m_state_new
    // when use_phase_frame = 0)
    const GRLevelData &phys_state = physical_state();
    
    Potential potential(m_p.potential_params);
    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                           phys_state, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                       phys_state, m_state_diagnostics, EXCLUDE_GHOST_CELLS);


#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        // Hack: if avg radius of found AH is negative, we reset initial guess
        double current_avg_AH_radius =  m_bh_amr.m_ah_finder.get(0)->get_ave_F();
        if (current_avg_AH_radius < 0.)
        {
            m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
            pout() << "AHFinder: resetting initial guess as avg radius is "
                      "negative."
                   << endl;
        }
        else if (current_avg_AH_radius > 10.)
        {
            m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
            pout() << "AHFinder: resetting initial guess as avg radius is "
                      "too large."
                   << endl;
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif

    if (m_p.activate_extraction == 1 &&
       at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
            {
                pout() << "BinaryBSLevel::specificPostTimeStep:"
                          " Extracting gravitational waves." << endl;
            }


            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(NoetherCharge<FourthOrderDerivatives>(m_dx, m_p.axis_corrected_charge), phys_state, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
        BoxLoops::loop(ADMQuantities(m_p.extraction_params.center, m_dx, 0), phys_state, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
	BoxLoops::loop(Constraints<Potential>(m_dx, potential, m_p.m_G_Newton), phys_state, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
    }

    if (m_level == 0)
    {
        bool first_step = (m_time == 0.);
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom2), 2, true);
        SmallDataIO constraints_file("constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});

        int adm_min_level = 0;
        bool calculate_adm = at_level_timestep_multiple(adm_min_level);
        if (calculate_adm)
        {   
            AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
            double M_ADM = amr_reductions.sum(c_Madm);
            SmallDataIO M_ADM_file("M_ADM", m_dt, m_time,
                                m_restart_time, SmallDataIO::APPEND,
                                first_step);
            M_ADM_file.remove_duplicate_time_data();
            if (first_step)
            {
                M_ADM_file.write_header_line({"M_ADM"});
            }
            M_ADM_file.write_time_data_line({M_ADM});
        }

        double noether_charge = amr_reductions.sum(c_N);
        SmallDataIO noether_charge_file("NoetherCharge", m_dt, m_time,
                                        m_restart_time,
                                        SmallDataIO::APPEND,
                                        first_step);
        noether_charge_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            noether_charge_file.write_header_line({"Noether Charge"});
        }
        noether_charge_file.write_time_data_line({noether_charge});

        AMRReductions<VariableType::evolution> amr_reductions_evolution(m_gr_amr);
	double min_chi = amr_reductions_evolution.min(c_chi);
        SmallDataIO min_chi_file("min_chi",
                                        m_dt, m_time, m_restart_time,
                                        SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (first_step)
        {
            min_chi_file.write_header_line({"min_chi"});
        }
        min_chi_file.write_time_data_line({min_chi});

        // Compute the maximum of mod_phi and write it to a file
        double mod_phi_max = amr_reductions.max(c_mod_phi);
        SmallDataIO mod_phi_max_file("mod_phi_max", m_dt, m_time,
                                 m_restart_time,
                                 SmallDataIO::APPEND,
                                 first_step);
        mod_phi_max_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            mod_phi_max_file.write_header_line({"max mod phi"});
        }
        mod_phi_max_file.write_time_data_line({mod_phi_max});

        // Phase frame: range of the frame frequency and phase over all levels
        if (m_p.phase_frame_params.use_phase_frame &&
            m_p.phase_frame_params.write_diagnostics)
        {
            const double omega_min = amr_reductions_evolution.min(c_frame_omega);
            const double omega_max = amr_reductions_evolution.max(c_frame_omega);
            const double theta_min = amr_reductions_evolution.min(c_frame_theta);
            const double theta_max = amr_reductions_evolution.max(c_frame_theta);
            SmallDataIO frame_file("frame_diagnostics", m_dt, m_time,
                                   m_restart_time, SmallDataIO::APPEND,
                                   first_step);
            frame_file.remove_duplicate_time_data();
            if (first_step)
            {
                frame_file.write_header_line(
                    {"min Omega", "max Omega", "min Theta", "max Theta"});
            }
            frame_file.write_time_data_line(
                {omega_min, omega_max, theta_min, theta_max});
        }

    }
    
}

const GRLevelData &ScalarField2DLevel::physical_state()
{
    if (!m_p.phase_frame_params.use_phase_frame)
        return m_state_new;
    frame_to_physical(m_state_new, m_state_phys);
    return m_state_phys;
}

void ScalarField2DLevel::frame_to_physical(const GRLevelData &a_frame,
                                           GRLevelData &a_phys) const
{
    a_phys.define(a_frame.disjointBoxLayout(), a_frame.nComp(),
                  a_frame.ghostVect());
    // copy everything, including ghost cells
    for (DataIterator dit = a_frame.dataIterator(); dit.ok(); ++dit)
    {
        a_phys[dit()].copy(a_frame[dit()]);
    }
    // rotate the matter slots, phi = exp(i Theta) psi, pointwise
    BoxLoops::loop(PhaseFrameRotate(+1), a_phys, a_phys, INCLUDE_GHOST_CELLS);
}
