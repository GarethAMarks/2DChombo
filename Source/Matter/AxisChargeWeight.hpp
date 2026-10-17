/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AXISCHARGEWEIGHT_HPP_
#define AXISCHARGEWEIGHT_HPP_

#include <cmath>

// Axis correction of cartoon (axisymmetric) cell sums  sum_j 2 pi y_j f_j dx^2.
//
// On the cell-centred grid y_j = (j + 1/2) dx, the sum is the midpoint rule for
// int_0^inf g(y) dy with g = y F(y), F even in y. Since g'(0) = F(0) != 0, the
// Euler-Maclaurin formula gives
//   int_0^inf g dy = dx sum_j g(y_j) - (dx^2 / 24) g'(0) + O(dx^4),
// i.e. the plain sum is only 2nd-order accurate, and a conserved integral
// (such as the Noether charge) shows a spurious O(dx^2) drift whenever the
// integrand changes near the axis. Two corrected weight sets, both 4th-order
// accurate, are available (axis_weight_factor returns w_j / y_j):
//
//  AxisWeight::euler_maclaurin (1): F(0) = sum_{j<3} e_j F(y_j) + O(dx^6)
//    (even extrapolation, e = (75, -25/2, 3/2)/64), so
//      w_j = y_j - (dx / 24) e_j   (j = 0, 1, 2),   w_j = y_j otherwise.
//  AxisWeight::flux_form (2): F(0) = F(y_0) + O(dx^2) is enough at 4th order:
//      w_0 = y_0 - dx / 24 = (11/12) y_0,   w_j = y_j otherwise.
//    These are the weights of the charge-conserving flux-form Klein-Gordon
//    operator (FluxFormKGCartoon.hpp): with them the operator is exactly
//    self-adjoint and 4th-order accurate up to the axis, so they must be used
//    for the charge (and the KO charge projection) whenever KG_flux_form = 1.
namespace AxisWeight
{
enum
{
    none = 0,
    euler_maclaurin = 1,
    flux_form = 2
};
}

inline double axis_weight_factor(const double y, const double dx,
                                 const int mode)
{
    static const double e[3] = {1.171875, -0.1953125, 0.0234375};
    if (mode == AxisWeight::none)
        return 1.0;
    const int j = static_cast<int>(std::floor(y / dx + 1e-9));
    if (mode == AxisWeight::flux_form)
        return (j == 0) ? 11.0 / 12.0 : 1.0;
    if (j < 0 || j > 2)
        return 1.0;
    return 1.0 - dx * e[j] / (24.0 * y);
}

/// Euler-Maclaurin weights (backwards-compatible two-argument form)
inline double axis_weight_factor(const double y, const double dx)
{
    return axis_weight_factor(y, dx, AxisWeight::euler_maclaurin);
}

#endif /* AXISCHARGEWEIGHT_HPP_ */
