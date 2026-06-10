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

// Star tracking
#include "STAMR.hpp"

#ifdef USE_TWOPUNCTURES
#include "TPAMR.hpp"
TPAMR tp_amr;
#endif


int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    #ifdef USE_TWOPUNCTURES
        //TPAMR bh_amr;
        tp_amr.set_two_punctures_parameters(sim_params.tp_params);
        // Run TwoPunctures solver if id_choice is appropriate
        if (sim_params.bosonstar_params.id_choice > 0)
            tp_amr.m_two_punctures.Run();
    #endif

    //BHAMR bh_amr;




    STAMR st_amr;

    // must be before 'setupAMRObject' to define punctures for tagging criteria
    if (sim_params.do_puncture_track)
    {
        // the tagging criterion used in this example means that the punctures
        // should be on the max level but let's fill ghosts on the level below
        // too just in case
        int puncture_tracker_min_level = sim_params.max_level - 1;
        st_amr.m_puncture_tracker.initial_setup(
            {sim_params.bh1_params.center, sim_params.bh2_params.center},
            "punctures", sim_params.data_path, puncture_tracker_min_level);
    }

    if (sim_params.do_star_track)
    {
        st_amr.m_star_tracker.initialise_star_tracking(
            sim_params.number_of_stars,
            {sim_params.positionA, sim_params.positionB},
            sim_params.star_points, sim_params.star_track_width_A,
            sim_params.star_track_width_B,
            sim_params.star_track_direction_of_motion);
    }

    DefaultLevelFactory<ScalarField2DLevel> scalarfield2D_level_fact(st_amr, sim_params);
    setupAMRObject(st_amr, scalarfield2D_level_fact);

    // call this after amr object setup so grids known
    // and need it to stay in scope throughout run
    AMRInterpolator<Lagrange<4>> interpolator(
        st_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    st_amr.set_interpolator(&interpolator);


    // must be after interpolator is set
    if (sim_params.do_star_track)
        st_amr.m_star_tracker.restart_star_tracking();

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
    call_task.execute(st_amr);

    st_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    st_amr.conclude();

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
