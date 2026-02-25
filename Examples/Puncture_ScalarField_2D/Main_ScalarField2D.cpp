/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CH_Timer.H"
#include "parstream.H" //Gives us pout()
#include <chrono>
#include <iostream>

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "ScalarField2DLevel.hpp"

#ifdef USE_TWOPUNCTURES
#include "TPAMR.hpp"
//TPAMR tp_amr;
#endif


int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    pout() << "starting param parsing" << endl;
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

   // pout() << "finished param parsing" << endl;
    


    BHAMR bh_amr;

        // must be before 'setupAMRObject' to define punctures for tagging criteria
    if (sim_params.do_puncture_track)
    {
        // the tagging criterion used in this example means that the punctures
        // should be on the max level but let's fill ghosts on the level below
        // too just in case
        int puncture_tracker_min_level = sim_params.max_level - 1;
        bh_amr.m_puncture_tracker.initial_setup(
            {sim_params.bh1_params.center, sim_params.bh2_params.center},
            "punctures", sim_params.data_path, puncture_tracker_min_level,
            sim_params.coarsest_dt);
                // Compute initial star positions in simulation coordinates.
        // BosonStar.impl.hpp measures coords relative to star_centre
        // (the centre of mass), with:
        //   Star 1 (rapidity > 0, moving right): offset +q*d/(q+1) in x, -q*b/(q+1) in y
        //   Star 2 (rapidity2, moving left):     offset  -d/(q+1)  in x,   +b/(q+1) in y
        // where d = BS_separation, b = BS_impact_parameter, q = mass_ratio.

        /*const auto &bs1 = sim_params.bosonstar_params;
        const double q = bs1.mass_ratio;
        const double d = bs1.BS_separation;
        const double b = bs1.BS_impact_parameter;
        const auto &com = bs1.star_centre; // centre of mass in sim coords

        std::array<double, CH_SPACEDIM> star1_pos = com;
        star1_pos[0] += q * d / (q + 1.);
        star1_pos[1] -= q * b / (q + 1.);

        std::array<double, CH_SPACEDIM> star2_pos = com;
        star2_pos[0] -= d / (q + 1.);
        star2_pos[1] += b / (q + 1.);

        bh_amr.m_puncture_tracker.initial_setup(
            {star1_pos, star2_pos},
            "punctures", sim_params.data_path, puncture_tracker_min_level);*/
    }


    DefaultLevelFactory<ScalarField2DLevel> scalarfield2D_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, scalarfield2D_level_fact);

    // call this after amr object setup so grids known
    // and need it to stay in scope throughout run
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    
    bh_amr.set_interpolator(&interpolator);

    // Initialise puncture tracker: writes t=0 header+data to punctures.dat
    // and sets m_num_punctures. Without this, execute_tracking() is a no-op.
    if (sim_params.do_puncture_track)
    {
        bh_amr.m_puncture_tracker.restart_punctures();
    }


#ifdef USE_AHFINDER
    if (sim_params.AH_activate)
    {
        AHSurfaceGeometry sph1(sim_params.bh1_params.center);
        // AHSurfaceGeometry sph2(sim_params.bh2_params.center);

        bh_amr.m_ah_finder.add_ah(sph1, sim_params.AH_1_initial_guess_ellipsoid,
                                  sim_params.AH_params);
        // bh_amr.m_ah_finder.add_ah(sph2, sim_params.AH_2_initial_guess_ellipsoid,
        //                           sim_params.AH_params);
        // bh_amr.m_ah_finder.add_ah_merger(0, 1, sim_params.AH_params);
    }
#endif

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    // Add a scheduler to call specificPostTimeStep on every AMRLevel at t=0
    auto task = [](GRAMRLevel *level) {
        if (level->time() == 0.)
            level->specificPostTimeStep();
    };
    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(bh_amr);

    bh_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bh_amr.conclude();

    CH_TIMER_REPORT(); // Report results when running with Chombo timers.

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
