/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4CartoonUserVariables.hpp"
#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum

    // Phase frame (see Source/PhaseFrame): the field is stored as
    // psi = exp(-i frame_theta) phi in the phi/phi_Im slots (and likewise
    // P = exp(-i frame_theta) Pi in the Pi/Pi_Im slots), with
    // d_t frame_theta = frame_omega. With use_phase_frame = 0 both stay zero
    // and the slots hold the physical field as before. (The name "Theta" is
    // taken by the CCZ4 constraint variable.)
    c_frame_theta = NUM_CCZ4_VARS,
    c_frame_omega,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"frame_theta", "frame_omega"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
