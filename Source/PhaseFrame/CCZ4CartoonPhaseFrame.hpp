/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4CARTOONPHASEFRAME_HPP_
#define CCZ4CARTOONPHASEFRAME_HPP_

#include "CCZ4Cartoon.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "PhaseFrameParams.hpp"
#include "PhaseFrameVars.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"
#include <cmath>

/// CCZ4 cartoon right-hand side with the scalar field in a U(1) phase frame
/**
 * Evolved matter variables: psi = exp(-i Theta) phi, P = exp(-i Theta) Pi
 * (stored in the phi/Pi slots), plus the frame phase Theta (frame_theta) and
 * frequency Omega (frame_omega):
 *
 *   d_t psi   = exp(-i Theta) [d_t phi] - i Omega psi       (same for P)
 *   d_t Theta = Omega
 *   d_t Omega = eta Im(conj(psi) d_t psi) / (|psi|^2 + eps^2)
 *               + kappa lap Omega + zeta lap Theta          (local driver)
 *   d_t Omega = 0                                           (fixed frame)
 *
 * where [d_t phi] is the unchanged CCZ4Cartoon scalar-field right-hand side
 * evaluated in the local frame, and lap is the flat cartoon Laplacian
 * d_xx + d_yy + nS d_y / y (nS = GR_SPACEDIM - CH_SPACEDIM). kappa =
 * kappa_fac * dx on each level. Theta and Omega get no KO dissipation.
 *
 * Theta(x, y) is an axisymmetric scalar and the frame is a rotation in field
 * space, so it commutes with the modified cartoon reduction: the cartoon terms
 * (e.g. -lapse chi h^ww d_y phi / y in the Pi equation) are evaluated on the
 * physical field exactly as before. Two ways of evaluating [d_t phi]:
 *
 *  links = true:  the input box must hold the physical field
 *                 exp(i Theta) psi in the matter slots (and the unchanged
 *                 Theta, Omega), see ScalarField2DLevel::specificEvalRHS.
 *                 All stencils (derivatives, advection, KO) then act on the
 *                 physical field, i.e. on exp(i(Theta_k - Theta_j)) psi_k in
 *                 the frame of the central cell j: the semi-discretisation is
 *                 identical to the original one and only the variables
 *                 integrated in time (and interpolated between levels) change.
 *  links = false: the input box holds psi; derivatives are replaced by the
 *                 continuum covariant ones,
 *                   D_i psi    = d_i psi + i (d_i Theta) psi,
 *                   D_i D_j psi = d_ij psi + i (Theta_i d_j psi + Theta_j d_i
 *                                 psi) + i Theta_ij psi - Theta_i Theta_j psi,
 *                   shift^i D_i psi = advec psi + i (shift.d Theta) psi,
 *                 and KO dissipation acts on psi directly.
 * The stress tensor is U(1) invariant and is the same in either case.
 */
template <class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives,
          class potential_t = Potential>
class CCZ4CartoonPhaseFrame
    : public CCZ4Cartoon<gauge_t, deriv_t, potential_t>
{
  public:
    using BaseCartoon = CCZ4Cartoon<gauge_t, deriv_t, potential_t>;
    using params_t = typename BaseCartoon::params_t;

    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsWithGauge<data_t>;

    CCZ4CartoonPhaseFrame(params_t a_params, double a_dx, double a_sigma,
                          potential_t a_potential, double a_G_Newton,
                          const PhaseFrame_params_t &a_frame_params,
                          int a_formulation,
                          double a_cosmological_constant = 0.)
        : BaseCartoon(a_params, a_dx, a_sigma, a_potential, a_G_Newton,
                      a_formulation, a_cosmological_constant),
          m_frame_params(a_frame_params),
          m_kappa(a_frame_params.kappa_fac * a_dx)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        using namespace PhaseFrame;

        // Frame variables are never rotated: same in either input box
        const auto fvars = current_cell.template load_vars<FrameVars>();
        const auto d1f = this->m_deriv.template diff1<FrameVars>(current_cell);

        auto vars = current_cell.template load_vars<Vars>();
        auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
        auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
        auto advec =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

        Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

        const data_t cos_th = cos(fvars.theta);
        const data_t sin_th = sin(fvars.theta);

        if (!m_frame_params.links)
        {
            // continuum connection: needs second derivatives of Theta and the
            // (upwinded) advection of Theta
            const auto d2f =
                this->m_deriv.template diff2<FrameVars>(current_cell);
            const auto advec_f = this->m_deriv.template advection<FrameVars>(
                current_cell, vars.shift);
            apply_connection(vars, d1, d2, advec, d1f.theta, d2f.theta,
                             advec_f.theta);
        }

        // Unchanged cartoon CCZ4 + scalar field right-hand side (in the
        // local frame), and KO on whatever the input box holds
        Vars<data_t> rhs;
        this->rhs_equation(rhs, vars, d1, d2, advec, coords.y);
        // flux-form Klein-Gordon operator (KG_flux_form; links mode only, on
        // the physical field)
        this->add_KG_flux_form(rhs, current_cell, coords.y);
        this->m_deriv.add_dissipation(rhs, current_cell, this->m_sigma);

        // psi, P and the local-frame time derivatives exp(-i Theta) d_t phi
        data_t psi, psi_Im, P, P_Im, dpsi, dpsi_Im, dP, dP_Im;
        if (m_frame_params.links)
        {
            // input holds the physical field: rotate by exp(-i Theta)
            const data_t s = -sin_th;
            complex_rotate(vars.phi, vars.phi_Im, cos_th, s, psi, psi_Im);
            complex_rotate(vars.Pi, vars.Pi_Im, cos_th, s, P, P_Im);
            complex_rotate(rhs.phi, rhs.phi_Im, cos_th, s, dpsi, dpsi_Im);
            complex_rotate(rhs.Pi, rhs.Pi_Im, cos_th, s, dP, dP_Im);
        }
        else
        {
            psi = vars.phi;
            psi_Im = vars.phi_Im;
            P = vars.Pi;
            P_Im = vars.Pi_Im;
            dpsi = rhs.phi;
            dpsi_Im = rhs.phi_Im;
            dP = rhs.Pi;
            dP_Im = rhs.Pi_Im;
        }

        // rotating-frame term -i Omega psi
        const data_t &omega = fvars.omega;
        dpsi += omega * psi_Im;
        dpsi_Im -= omega * psi;
        dP += omega * P_Im;
        dP_Im -= omega * P;

        rhs.phi = dpsi;
        rhs.phi_Im = dpsi_Im;
        rhs.Pi = dP;
        rhs.Pi_Im = dP_Im;

        // Frame evolution
        data_t rhs_omega = 0.;
        if (m_frame_params.local_driver)
        {
            const data_t mod_psi2 = psi * psi + psi_Im * psi_Im;
            const data_t im_psibar_dpsi = psi * dpsi_Im - psi_Im * dpsi;
            const double eps2 = m_frame_params.eps * m_frame_params.eps;
            rhs_omega = m_frame_params.eta * im_psibar_dpsi / (mod_psi2 + eps2);

            if (m_kappa != 0. || m_frame_params.stiffness != 0.)
            {
                // second derivatives of the frame variables (links mode has
                // not computed them yet)
                const auto d2f_lap =
                    this->m_deriv.template diff2<FrameVars>(current_cell);
                if (m_kappa != 0.)
                {
                    rhs_omega += m_kappa * cartoon_laplacian(d1f.omega,
                                                             d2f_lap.omega,
                                                             coords.y);
                }
                if (m_frame_params.stiffness != 0.)
                {
                    rhs_omega += m_frame_params.stiffness *
                                 cartoon_laplacian(d1f.theta, d2f_lap.theta,
                                                   coords.y);
                }
            }
        }

        current_cell.store_vars(rhs);
        current_cell.store_vars(omega, c_frame_theta);
        current_cell.store_vars(rhs_omega, c_frame_omega);
    }

  protected:
    const PhaseFrame_params_t m_frame_params;
    const double m_kappa; //!< Omega smoothing coefficient on this level

    /// Replace the matter derivatives of psi by the U(1)-covariant ones
    /// (continuum connection, links = false)
    template <class data_t>
    static void apply_connection(const Vars<data_t> &vars,
                                 Vars<Tensor<1, data_t>> &d1,
                                 Diff2Vars<Tensor<2, data_t>> &d2,
                                 Vars<data_t> &advec,
                                 const Tensor<1, data_t> &dtheta,
                                 const Tensor<2, data_t> &d2theta,
                                 const data_t &advec_theta)
    {
        // second derivatives first: they use the uncorrected d1
        FOR(i, j)
        {
            const data_t a_ij = d2.phi[i][j];
            const data_t b_ij = d2.phi_Im[i][j];
            d2.phi[i][j] = a_ij - dtheta[i] * d1.phi_Im[j] -
                           dtheta[j] * d1.phi_Im[i] -
                           d2theta[i][j] * vars.phi_Im -
                           dtheta[i] * dtheta[j] * vars.phi;
            d2.phi_Im[i][j] = b_ij + dtheta[i] * d1.phi[j] +
                              dtheta[j] * d1.phi[i] +
                              d2theta[i][j] * vars.phi -
                              dtheta[i] * dtheta[j] * vars.phi_Im;
        }
        FOR(i)
        {
            d1.phi[i] -= dtheta[i] * vars.phi_Im;
            d1.phi_Im[i] += dtheta[i] * vars.phi;
            d1.Pi[i] -= dtheta[i] * vars.Pi_Im;
            d1.Pi_Im[i] += dtheta[i] * vars.Pi;
        }
        advec.phi -= advec_theta * vars.phi_Im;
        advec.phi_Im += advec_theta * vars.phi;
        advec.Pi -= advec_theta * vars.Pi_Im;
        advec.Pi_Im += advec_theta * vars.Pi;
    }
};

#endif /* CCZ4CARTOONPHASEFRAME_HPP_ */
