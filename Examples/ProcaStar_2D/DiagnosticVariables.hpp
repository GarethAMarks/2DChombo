/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom, // sqrt(Mom_1^2 + Mom_2^2)

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_Madm,
    c_Padm,

    c_rho,
    c_rho_ADM,
    c_Sx,
    c_Sy,

    c_N,        // particle-number integrand n chi^{-3/2} 2 pi y
    c_mod_Xphi, // |X_phi|: equals A_central = f(0)/alpha(0) at the centre of a star at rest
    c_n,        // particle number density measured by the normal observers, Im(E^i Xbar_i) (n = -n_a j^a)
    c_Gauss,    // |D_i E^i + mu^2 X_phi|

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham", "Mom1", "Mom2", "Mom", "Weyl4_Re", "Weyl4_Im", "M_adm", "P_adm",
    "rho", "rho_ADM", "Sx", "Sy", "N", "mod_Xphi", "n", "Gauss"};
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
