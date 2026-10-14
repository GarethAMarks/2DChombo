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

    // sqrt(Mom_1^2 + Mom_2^2)
    c_Mom,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_Madm,
    c_Padm,

    c_rho,
    c_rho_ADM, // basically rho * sqrt(gamma)

    c_Sx,
    c_Sy,

    c_N,
    c_mod_phi,
    c_dt_mod_phi,
    c_gamma_tt,

    c_profile1,
    c_profile2,

    // physical (lab-frame) scalar field exp(i frame_theta) psi; equal to
    // phi/phi_Im when use_phase_frame = 0
    c_phi_phys,
    c_phi_Im_phys,

    // c_Sxx,
    // c_Sxy,
    // c_Syy,

    // c_Sww,

    // c_S,

    // c_Sxx_TF,
    // c_Sxy_TF,
    // c_Syy_TF,
    // c_Sww_TF,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
// NB: one name per enum entry above, in the same order ("dt_mod_phi" and
// "gamma_tt" were previously missing, which shifted "profile1"/"profile2" onto
// the wrong variables)
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",     "Mom",

    "Weyl4_Re", "Weyl4_Im",

    "M_adm",    "P_adm",

    "rho",      "rho_ADM",  "Sx",  "Sy",

    "N", "mod_phi", "dt_mod_phi", "gamma_tt",

    "profile1", "profile2",

    "phi_phys", "phi_Im_phys"

    };
// "Sxx",    "Sxy",    "Syy",
// "Sww",      "S",        "Sxx_TF", "Sxy_TF", "Syy_TF", "Sww_TF"
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
