/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BLACKSTRINGLEVEL_HPP_
#define BLACKSTRINGLEVEL_HPP_

#undef USE_AHFINDER //just ignore ah finder for now

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"

class ScalarField2DLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ScalarField2DLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

     // Typedef for scalar field
    typedef ComplexScalarField<Potential> ComplexScalarFieldWithPotential;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state,
                                         const double a_maximum) override;

    virtual void specificPostTimeStep() override;

    // U(1) phase frame helpers (see Source/PhaseFrame)

    /// Copy of the state with the scalar field rotated back to the physical
    /// frame, phi = exp(i frame_theta) psi (filled by physical_state())
    GRLevelData m_state_phys;

    /// m_state_new with the physical scalar field: m_state_new itself if
    /// use_phase_frame = 0, otherwise a rotated copy in m_state_phys. Call
    /// after fillAllGhosts(); used by every diagnostic that is not U(1)
    /// invariant (stress tensor, constraints, Weyl, tagging, plots)
    const GRLevelData &physical_state();

    /// a_phys = a_frame with matter rotated by exp(+i frame_theta), including
    /// ghost cells (a_phys is (re)defined on a_frame's layout)
    void frame_to_physical(const GRLevelData &a_frame,
                           GRLevelData &a_phys) const;

    /// Charge projection of the KO dissipation of the scalar field (see
    /// Source/Matter/KOChargeProjection.hpp): a_rhs -= lambda g(a_soln), with
    /// lambda such that the dissipation does not change the Noether charge on
    /// the cells of this level not covered by a finer level. a_rhs_input is the
    /// state the RHS and its KO term were computed from.
    void ko_charge_projection(const GRLevelData &a_rhs_input,
                              const GRLevelData &a_soln, GRLevelData &a_rhs);
};

#endif /* BLACKSTRINGLEVEL_HPP_ */
