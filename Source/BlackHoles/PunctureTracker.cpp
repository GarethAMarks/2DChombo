/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "PunctureTracker.hpp"
#include "GRAMR.hpp" // for static_cast to access get_restart_step()
#include "ChomboParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "InterpolationQuery.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data

//! Set punctures post restart
void PunctureTracker::initial_setup(
    const std::vector<std::array<double, CH_SPACEDIM>> &initial_puncture_coords,
    const std::string &a_filename, const std::string &a_output_path,
    const int a_min_level, const double a_coarsest_dt)
{
    if (!FilesystemTools::directory_exists(a_output_path))
        FilesystemTools::mkdir_recursive(a_output_path);

    m_punctures_filename = a_output_path + a_filename;

    // first set the puncture data
    // m_num_punctures is only set later
    m_puncture_coords = initial_puncture_coords;

    m_min_level = a_min_level;
    m_coarsest_dt = a_coarsest_dt;
}

void PunctureTracker::restart_punctures()
{
    // Use m_restart_step (set from checkpoint) rather than s_step
    // (which is the live global counter and may be ahead of the
    // restart step on a re-restart).
    const AMR &amr = m_interpolator->getAMR();
    // m_restart_step == -1 means fresh start (no checkpoint read)
    int restart_step = static_cast<const GRAMR &>(amr).get_restart_step();

    if (restart_step <= 0)
    {
        // Fresh start: use the param values supplied to initial_setup
        set_initial_punctures();
    }
    else
    {
        // Restarted run: read puncture location from the output file
        read_in_punctures(restart_step,
                          amr.getCurrentTime(), m_coarsest_dt);
    }
}

//! set and write initial puncture locations
void PunctureTracker::set_initial_punctures()
{
    CH_assert(m_puncture_coords.size() > 0); // sanity check

    m_num_punctures = m_puncture_coords.size();
    m_puncture_shift.resize(m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        // assume initial shift is always zero
        FOR(i) { m_puncture_shift[ipuncture][i] = 0.0; }
    }

    // now the write out to a new file
    bool first_step = true;
    double dt = 1.; // doesn't matter
    double time = 0.;
    double restart_time = 0.;
    SmallDataIO punctures_file(m_punctures_filename, dt, time, restart_time,
                               SmallDataIO::APPEND, first_step);
    std::vector<std::string> header1_strings(CH_SPACEDIM * m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        std::string idx = std::to_string(ipuncture + 1);
        header1_strings[CH_SPACEDIM * ipuncture + 0] = "x_" + idx;
        header1_strings[CH_SPACEDIM * ipuncture + 1] = "y_" + idx;
#if CH_SPACEDIM == 3
        header1_strings[CH_SPACEDIM * ipuncture + 2] = "z_" + idx;
#endif
    }
    punctures_file.write_header_line(header1_strings);

    // use a vector for the write out
    punctures_file.write_time_data_line(get_puncture_vector());
}

//! Set punctures post restart
void PunctureTracker::read_in_punctures(int a_int_step, double a_current_time,
                                        double a_coarsest_dt)
{
    // read them in from the Punctures file at current time m_time
    // NB opening in APPEND mode allows reading where m_restart_time
    // is greater than zero and m_time < m_restart_time + m_dt.
    // We must use the coarsest-level dt here so that SmallDataIO
    // opens the file in read+append mode (needed for duplicate removal),
    // and so that the coords_epsilon matching is correct.
    bool first_step = false;
    double dt = (a_coarsest_dt > 0.) ? a_coarsest_dt
                                     : (a_current_time / a_int_step);
    SmallDataIO punctures_file(m_punctures_filename, dt, a_current_time,
                               a_current_time, SmallDataIO::APPEND, first_step);

    // NB need to give the get function an empty vector to fill
    std::vector<double> puncture_vector;
    punctures_file.get_specific_data_line(puncture_vector, a_current_time);

    // check the data returned is the right size
    CH_assert(puncture_vector.size() % CH_SPACEDIM == 0);

    m_num_punctures = puncture_vector.size() / CH_SPACEDIM;
    m_puncture_coords.resize(m_num_punctures);

    // remove any duplicate data from the file
    const bool keep_m_time_data = true;
    punctures_file.remove_duplicate_time_data(keep_m_time_data);

    // convert vector to list of coords
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
#if CH_SPACEDIM == 2
        m_puncture_coords[ipuncture] = {
            puncture_vector[ipuncture * CH_SPACEDIM + 0],
            puncture_vector[ipuncture * CH_SPACEDIM + 1]};
#endif

#if CH_SPACEDIM == 3
        m_puncture_coords[ipuncture] = {
            puncture_vector[ipuncture * CH_SPACEDIM + 0],
            puncture_vector[ipuncture * CH_SPACEDIM + 1],
            puncture_vector[ipuncture * CH_SPACEDIM + 2]};
#endif
    }

    // set the coordinates and get the current shift
    interp_shift();

    // print out values into pout files
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        pout() << "Puncture " << ipuncture
               << " restarted at : " << m_puncture_coords[ipuncture][0] << " "
               << m_puncture_coords[ipuncture][1] << " "
               << m_puncture_coords[ipuncture][2] << endl;
        pout() << " with shift vector : " << m_puncture_shift[ipuncture][0]
               << " " << m_puncture_shift[ipuncture][1] << " "
               << m_puncture_shift[ipuncture][2] << endl;
        pout() << "at time = " << a_current_time << endl;
    }
}

//! Execute the tracking and write out
void PunctureTracker::execute_tracking(double a_time, double a_restart_time,
                                       double a_dt, const bool write_punctures)
{
    CH_TIME("PunctureTracker::execute_tracking");
    // leave if this is called at t=0, we don't want to move the puncture yet
    if (m_num_punctures == 0 || a_time == 0.)
        return;
    CH_assert(m_interpolator != nullptr); // sanity check

    // get puncture coordinates and old shift value
    std::vector<std::array<double, CH_SPACEDIM>> old_shift = m_puncture_shift;
    CH_assert(m_puncture_coords.size() == m_num_punctures); // sanity check

    // new shift value
    interp_shift();

    // update puncture locations using second order update
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        FOR(i)
        {
            m_puncture_coords[ipuncture][i] +=
                -0.5 * a_dt *
                (m_puncture_shift[ipuncture][i] + old_shift[ipuncture][i]);
        }
    }

    // print them out
    if (write_punctures)
    {
        bool first_step = false;
        // Use m_coarsest_dt (not the fine-level a_dt) so that the
        // SmallDataIO restart-window condition
        //   m_time < m_restart_time + m_dt + epsilon
        // fires correctly and remove_duplicate_time_data works.
        double dt_for_io = (m_coarsest_dt > 0.) ? m_coarsest_dt : a_dt;
        SmallDataIO punctures_file(m_punctures_filename, dt_for_io, a_time,
                                   a_restart_time, SmallDataIO::APPEND,
                                   first_step);
        punctures_file.remove_duplicate_time_data();

        // use a vector for the write out
        punctures_file.write_time_data_line(get_puncture_vector());
    }
}

//! Use the interpolator to get the value of the shift at
//! given coords
void PunctureTracker::interp_shift()
{
    CH_TIME("PunctureTracker::interp_shift");
    // resize the vector to the number of punctures
    m_puncture_shift.resize(m_num_punctures);

    // refresh interpolator
    bool fill_ghosts = false;
    m_interpolator->refresh(fill_ghosts);
    // only fill the ghosts we need
#if CH_SPACEDIM == 2
    m_interpolator->fill_multilevel_ghosts(
        VariableType::evolution, Interval(c_shift1, c_shift2), m_min_level);
#endif
#if CH_SPACEDIM == 3
    m_interpolator->fill_multilevel_ghosts(
        VariableType::evolution, Interval(c_shift1, c_shift3), m_min_level);
#endif

    // set up shift and coordinate holders
    std::vector<double> interp_shift1(m_num_punctures);
    std::vector<double> interp_shift2(m_num_punctures);
#if CH_SPACEDIM == 3
    std::vector<double> interp_shift3(m_num_punctures);
#endif
    std::vector<double> interp_x(m_num_punctures);
    std::vector<double> interp_y(m_num_punctures);
#if CH_SPACEDIM == 3
    std::vector<double> interp_z(m_num_punctures);
#endif

    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_x[ipuncture] = m_puncture_coords[ipuncture][0];
        interp_y[ipuncture] = m_puncture_coords[ipuncture][1];
#if CH_SPACEDIM == 3
        interp_z[ipuncture] = m_puncture_coords[ipuncture][2];
#endif
    }

    // setup query
    InterpolationQuery query(m_num_punctures);
#if CH_SPACEDIM == 2
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .addComp(c_shift1, interp_shift1.data())
        .addComp(c_shift2, interp_shift2.data());
#endif

#if CH_SPACEDIM == 3
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(c_shift1, interp_shift1.data())
        .addComp(c_shift2, interp_shift2.data())
        .addComp(c_shift3, interp_shift3.data());
#endif

    // engage!
    m_interpolator->interp(query);

    // put the shift values into the output array
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
#if CH_SPACEDIM == 2
        m_puncture_shift[ipuncture] = {interp_shift1[ipuncture],
                                       interp_shift2[ipuncture]};
#endif

#if CH_SPACEDIM == 3
        m_puncture_shift[ipuncture] = {interp_shift1[ipuncture],
                                       interp_shift2[ipuncture],
                                       interp_shift3[ipuncture]};
#endif
    }
}

//! get a vector of the puncture coords - used for write out
std::vector<double> PunctureTracker::get_puncture_vector() const
{
    std::vector<double> puncture_vector;
    puncture_vector.resize(m_num_punctures * CH_SPACEDIM);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        // We can probably replace this with a FOR(i)
        puncture_vector[ipuncture * CH_SPACEDIM + 0] =
            m_puncture_coords[ipuncture][0];
        puncture_vector[ipuncture * CH_SPACEDIM + 1] =
            m_puncture_coords[ipuncture][1];
#if CH_SPACEDIM == 3
        puncture_vector[ipuncture * CH_SPACEDIM + 2] =
            m_puncture_coords[ipuncture][2];
#endif
    }
    return puncture_vector;
}
