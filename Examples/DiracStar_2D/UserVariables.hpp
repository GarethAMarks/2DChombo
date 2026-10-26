/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

// This example evolves a Dirac field instead of the complex scalar field, so the scalar variables are left out of the
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
    // Dirac spinor (frame components on the cartoon plane, see DiracFieldCartoon.hpp):
    // psi0 = upper spin-up (even in y), psi1 = upper spin-down (odd), psi2 = lower spin-up (even), psi3 = lower
    // spin-down (odd); spin quantised along the symmetry axis x
    c_psi0_Re = NUM_CCZ4_VARS,
    c_psi0_Im,
    c_psi1_Re,
    c_psi1_Im,
    c_psi2_Re,
    c_psi2_Im,
    c_psi3_Re,
    c_psi3_Im,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS> user_variable_names = {
    "psi0_Re", "psi0_Im", "psi1_Re", "psi1_Im", "psi2_Re", "psi2_Im", "psi3_Re", "psi3_Im"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
