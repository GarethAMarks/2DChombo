/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHASEFRAMEVARS_HPP_
#define PHASEFRAMEVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp" // needs c_frame_theta, c_frame_omega
#include "VarsTools.hpp"

namespace PhaseFrame
{
/// Frame phase Theta and frame frequency Omega = d_t Theta
template <class data_t> struct FrameVars
{
    data_t theta;
    data_t omega;

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools;
        define_enum_mapping(mapping_function, c_frame_theta, theta);
        define_enum_mapping(mapping_function, c_frame_omega, omega);
    }
};

/// The four real matter components (rotated by the frame)
template <class data_t> struct MatterVars
{
    data_t phi;
    data_t phi_Im;
    data_t Pi;
    data_t Pi_Im;

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools;
        define_enum_mapping(mapping_function, c_phi, phi);
        define_enum_mapping(mapping_function, c_phi_Im, phi_Im);
        define_enum_mapping(mapping_function, c_Pi, Pi);
        define_enum_mapping(mapping_function, c_Pi_Im, Pi_Im);
    }
};

/// (re + i im) * (c + i s)
template <class data_t>
ALWAYS_INLINE void complex_rotate(const data_t &re, const data_t &im,
                                  const data_t &c, const data_t &s,
                                  data_t &out_re, data_t &out_im)
{
    out_re = re * c - im * s;
    out_im = re * s + im * c;
}

/// Flat Laplacian of an axisymmetric scalar in the cartoon reduction:
/// d_xx f + d_yy f + nS (d_y f)/y, y = cartoon (cylindrical) radius,
/// nS = GR_SPACEDIM - CH_SPACEDIM transverse dimensions. Regular on the axis
/// since d_y f is odd in y (and y >= dx/2 on the cell-centred grid).
template <class data_t>
ALWAYS_INLINE data_t cartoon_laplacian(const Tensor<1, data_t> &d1f,
                                       const Tensor<2, data_t> &d2f,
                                       const double cartoon_coord)
{
    const int dI = CH_SPACEDIM - 1;
    const int nS = GR_SPACEDIM - CH_SPACEDIM;
    data_t lap = nS * d1f[dI] / cartoon_coord;
    FOR(i) { lap += d2f[i][i]; }
    return lap;
}
} // namespace PhaseFrame

#endif /* PHASEFRAMEVARS_HPP_ */
