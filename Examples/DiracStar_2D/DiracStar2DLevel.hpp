/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIRACSTAR2DLEVEL_HPP_
#define DIRACSTAR2DLEVEL_HPP_

#undef USE_AHFINDER // just ignore ah finder for now

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

//! Level class for a Dirac star in axisymmetry (modified cartoon), the Dirac analogue of ScalarField2DLevel
class DiracStar2DLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<DiracStar2DLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);

    /// Things to do at every full timestep
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs, const double a_time) override;

    virtual void specificUpdateODE(GRLevelData &a_soln, const GRLevelData &a_rhs, Real a_dt) override;

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion, const FArrayBox &current_state,
                                         const double a_maximum) override;

    virtual void specificPostTimeStep() override;
};

#endif /* DIRACSTAR2DLEVEL_HPP_ */
