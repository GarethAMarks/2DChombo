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

    c_mode_power_re_0,
    c_mode_power_re_1,
    c_mode_power_re_2,
    c_mode_power_re_3,
    c_mode_power_re_4,
    c_mode_power_re_5,
    c_mode_power_re_6,
    c_mode_power_re_7,
    c_mode_power_re_8,
    c_mode_power_re_9,
    c_mode_power_re_10,
    c_mode_power_re_11,
    c_mode_power_re_12,

    c_mode_power_im_0,
    c_mode_power_im_1,
    c_mode_power_im_2,
    c_mode_power_im_3,
    c_mode_power_im_4,
    c_mode_power_im_5,
    c_mode_power_im_6,
    c_mode_power_im_7,
    c_mode_power_im_8,
    c_mode_power_im_9,
    c_mode_power_im_10,
    c_mode_power_im_11,
    c_mode_power_im_12,

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
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",     "Mom",

    "Weyl4_Re", "Weyl4_Im",

    "M_adm",    "P_adm",

    "rho",      "rho_ADM",  "Sx",  "Sy",

    "N", "mod_phi",

    "profile1", "profile2",

    "mode_power_re_0",  "mode_power_re_1",  "mode_power_re_2",
    "mode_power_re_3",  "mode_power_re_4",  "mode_power_re_5",
    "mode_power_re_6",  "mode_power_re_7",  "mode_power_re_8",
    "mode_power_re_9",  "mode_power_re_10", "mode_power_re_11",
    "mode_power_re_12",

    "mode_power_im_0",  "mode_power_im_1",  "mode_power_im_2",
    "mode_power_im_3",  "mode_power_im_4",  "mode_power_im_5",
    "mode_power_im_6",  "mode_power_im_7",  "mode_power_im_8",
    "mode_power_im_9",  "mode_power_im_10", "mode_power_im_11",
    "mode_power_im_12"
    
    };
// "Sxx",    "Sxy",    "Syy",
// "Sww",      "S",        "Sxx_TF", "Sxy_TF", "Syy_TF", "Sww_TF"
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
