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

    c_N,       // particle-number integrand n chi^{-3/2} 2 pi y (both spinors)
    c_mod_psi, // sqrt(4 pi psi^dag psi): equals A_central = g(0) at the centre of the initial star

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham", "Mom1", "Mom2", "Mom", "Weyl4_Re", "Weyl4_Im", "M_adm", "P_adm",
    "rho", "rho_ADM", "Sx", "Sy", "N", "mod_psi"};
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
