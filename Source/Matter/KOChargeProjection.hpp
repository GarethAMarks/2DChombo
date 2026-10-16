/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KOCHARGEPROJECTION_HPP_
#define KOCHARGEPROJECTION_HPP_

// Charge projection of the Kreiss-Oliger dissipation of the complex scalar
// field (cartoon, axisymmetric).
//
// KO acts on phi and Pi separately, which does not respect the Noether charge
//   N = sum over cells of  w q,   w = 2 pi y chi^{-3/2} dx^2,
//   q = phi_Im Pi - phi Pi_Im
// (the density integrated by NoetherCharge). KO therefore changes N at a
// constant rate even for a stationary star, and the change feeds back into the
// star's frequency (dephasing). The projection removes the charge-changing
// component of the dissipation: with D the KO contribution to the right-hand
// side and
//   g = (d phi, d phi_Im, d Pi, d Pi_Im) = (-Pi_Im, Pi, phi_Im, -phi)
//     ( = (i Pi, -i phi) in complex notation, the gradient of q, with
//       dq . g = |phi|^2 + |Pi|^2 > 0 ),
// the matter right-hand side gets  -lambda g  with
//   lambda = sum w (dq/du . D + q dw/dchi / w D_chi) / sum w (|phi|^2 + |Pi|^2),
// so that d N / dt from the dissipation (including KO of chi in the weight) is
// zero. lambda is a single number per level and RHS evaluation, computed over
// the cells of the level that are not covered by a finer level (the only cells
// that count in the composite charge). g is covariant under phase rotations,
// so the same correction applies in the U(1) phase frame: there the densities
// are computed from the field that KO actually acted on (the physical field
// with link variables, psi otherwise), and g from the evolved variables.
//
//  KOChargeRateDensity  - per cell: w (dq/du . D + q (dw/dchi)/w D_chi) and
//                         w (|phi|^2 + |Pi|^2) into components 0 and 1 of the
//                         output (constant factors 2 pi dx^2 dropped)
// The correction  rhs -= lambda g  itself is a pointwise update done in
// ScalarField2DLevel::ko_charge_projection.

#include "AxisChargeWeight.hpp"
#include "CCZ4CartoonVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"
#include <cmath>

namespace KOChargeProjection
{
/// The variables whose dissipation changes the charge: the matter fields and
/// the conformal factor in the weight
template <class data_t> struct Vars
{
    data_t phi;
    data_t phi_Im;
    data_t Pi;
    data_t Pi_Im;
    data_t chi;

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools;
        define_enum_mapping(mapping_function, c_phi, phi);
        define_enum_mapping(mapping_function, c_phi_Im, phi_Im);
        define_enum_mapping(mapping_function, c_Pi, Pi);
        define_enum_mapping(mapping_function, c_Pi_Im, Pi_Im);
        define_enum_mapping(mapping_function, c_chi, chi);
    }
};
} // namespace KOChargeProjection

/// Per-cell charge rate of the KO dissipation and normalisation density.
/// Loop over the box the RHS was evaluated on (the same input as CCZ4Cartoon),
/// writing into a 2-component output.
class KOChargeRateDensity
{
  public:
    KOChargeRateDensity(const double a_dx, const double a_sigma,
                        const bool a_include_chi,
                        const bool a_axis_correction = false)
        : m_deriv(a_dx), m_sigma(a_sigma), m_include_chi(a_include_chi),
          m_axis_correction(a_axis_correction)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto u = current_cell.template load_vars<KOChargeProjection::Vars>();

        // KO contribution exactly as CCZ4Cartoon adds it (same operator and sigma)
        KOChargeProjection::Vars<data_t> D;
        VarsTools::assign(D, 0.);
        m_deriv.add_dissipation(D, current_cell, m_sigma);

        Coordinates<data_t> coords(current_cell, m_deriv.m_dx);
        const data_t &chi = u.chi;                     // as in NoetherCharge
        data_t w = coords.y * pow(chi, -1.5);
        if (m_axis_correction) // same weights as the axis-corrected charge
            w = w * axis_weight_factor(coords.y, m_deriv.m_dx);

        const data_t q = u.phi_Im * u.Pi - u.phi * u.Pi_Im;
        // dq/du . D  for u = (phi, phi_Im, Pi, Pi_Im)
        data_t rate = D.phi_Im * u.Pi + u.phi_Im * D.Pi - D.phi * u.Pi_Im -
                      u.phi * D.Pi_Im;
        if (m_include_chi)
            rate -= 1.5 * q * D.chi / chi; // (dw/dchi) / w = -3 / (2 chi)
        rate *= w;

        const data_t norm = w * (u.phi * u.phi + u.phi_Im * u.phi_Im +
                                 u.Pi * u.Pi + u.Pi_Im * u.Pi_Im);

        current_cell.store_vars(rate, 0);
        current_cell.store_vars(norm, 1);
    }

  protected:
    const FourthOrderDerivatives m_deriv;
    const double m_sigma;
    const bool m_include_chi;
    const bool m_axis_correction;
};

#endif /* KOCHARGEPROJECTION_HPP_ */
