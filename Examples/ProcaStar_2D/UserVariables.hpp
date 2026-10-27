/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// This example evolves a Proca field instead of the complex scalar field, so the scalar variables are left out of the
// CCZ4 cartoon variables (see CCZ4CartoonUserVariables.hpp and ADMConformalVars.hpp)
#ifndef CARTOON_NO_SCALAR_FIELD
#define CARTOON_NO_SCALAR_FIELD
#endif

#include "ArrayTools.hpp"
#include "CCZ4CartoonUserVariables.hpp"
#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    // complex Proca field (see ProcaFieldCartoon.hpp): X_x, X_y, E^x, E^y, X_phi = -n^a A_a and the Gauss-constraint
    // damping variable Z (the azimuthal components vanish on the cartoon plane)
    c_Xx_Re = NUM_CCZ4_VARS,
    c_Xy_Re,
    c_Xx_Im,
    c_Xy_Im,
    c_Ex_Re,
    c_Ey_Re,
    c_Ex_Im,
    c_Ey_Im,
    c_Xphi_Re,
    c_Xphi_Im,
    c_Z_Re,
    c_Z_Im,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS> user_variable_names = {
    "Xx_Re", "Xy_Re", "Xx_Im", "Xy_Im", "Ex_Re", "Ey_Re", "Ex_Im", "Ey_Im", "Xphi_Re", "Xphi_Im", "Z_Re", "Z_Im"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
