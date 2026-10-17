/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLUXFORMKGCARTOON_HPP_
#define FLUXFORMKGCARTOON_HPP_

// Charge-conserving (flux-form) Klein-Gordon operator for the cartoon
// (axisymmetric, CH_SPACEDIM = 2) complex scalar field; used by CCZ4Cartoon
// when KG_flux_form = 1.
//
// The Pi equation contains the spatial operator acting on phi
//   L phi = -alpha chi h^ww d_y phi / y + chi alpha Gamma^i d_i phi
//           + h^ij (-chi d_j alpha d_i phi + 1/2 d_i chi alpha d_j phi
//                   - chi alpha d_i d_j phi),
// which (with the computed contracted Christoffels and det h = 1) is
//   L phi = -(1/w) d_a ( w alpha chi h^ab d_b phi ) = -(1/w) d_a (y B^ab d_b phi),
//   w = y chi^{-3/2}  (the Noether-charge weight),  B^ab = alpha chi^{-1/2} h^ab
// (a, b = x, y; h^ab the inverse of the 2D block). The product-rule
// discretisation makes W L symmetric only up to O(dx^4), so the discrete
// Noether charge  sum W (phi_Im Pi - phi Pi_Im)  is conserved only to O(dx^4)
// (a bounded oscillation for an oscillating star). Here W L = -M with M an
// exactly symmetric (and positive semi-definite) matrix on the half plane, so the
// spatial operator conserves the discrete charge to round-off:
//
//   M = Dx^T diag(yhat B^xx) Dx + Dy^T diag(|y| B^yy) Dy + C + C^T,
//   C = Dy^T diag(|y| B^yx) Iy dx,
//   W = yhat chi^{-3/2},   yhat_j = y_j  except  yhat_0 = (11/12) y_0,
//
// with Dx, Dy the 4th-order staggered derivatives (cells -> faces)
//   (D phi)_{k+1/2} = (phi_{k-1} - 27 phi_k + 27 phi_{k+1} - phi_{k+2}) / (24 dx),
// Iy the 4th-order interpolation cells -> y faces, dx the collocated 4th-order
// centred derivative and B interpolated to faces with the same 4-point formula.
// Written cell by cell:
//   W L phi = - [ stagdiv_x(yhat B^xx Dx phi) + stagdiv_y(|y_f| B^yy Dy phi)
//                 + stagdiv_y(|y_f| B^yx Iy dx phi) + dx Iy^T(|y_f| B^yx Dy phi) ]
// (Iy^T = the same 4-point interpolation from y faces to cell centres).
//
// Axis closure. The axis y = 0 is a cell face; the ghost cells below it are
// filled by the reflective boundary with the usual parities (h_xy odd, all other
// variables used here even). With |y_f| multiplying every y-face coefficient,
// the face on the axis carries no flux and the ghost-cell formulas above are
// exactly the transposes on the half plane y > 0, so M is symmetric. Accuracy
// at the axis comes from the weight of the first cell row: with the plain
// midpoint weight y_0 the first row has an O(1) error; with
// yhat_0 = (11/12) y_0 (which is also a 4th-order quadrature of the axis
// integral, see AxisChargeWeight.hpp) its error is O(dx^2) and the solution
// stays 4th-order accurate (Python prototype: errors 3.2e-5, 2.0e-6, 1.3e-7,
// 8.1e-9 for dx = 0.2 ... 0.025). The charge and the KO charge projection must
// use the same weights (AxisWeight::flux_form).
//
// Only alpha, chi and h_ab enter (no metric derivatives). The stencil reaches 3
// cells in each direction (corners (+-2, +-3)), within GRChombo's 3 ghosts. In
// the phase frame the operator must act on the physical field (frame_links = 1).
// Not included: hww (det h = 1 assumed, as in NoetherCharge's weight).

#include "AxisChargeWeight.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"
#include <cmath>

// (Written for CH_SPACEDIM = 2: the x and y directions of the cartoon plane.)

class FluxFormKGCartoon
{
  public:
    FluxFormKGCartoon(const double a_dx) : m_dx(a_dx) {}

    /// The flux-form operator L phi for phi (out_re) and phi_Im (out_im) at the
    /// current cell, to be added to rhs.Pi and rhs.Pi_Im in place of the
    /// product-rule terms. y = cartoon coordinate of the cell (> 0).
    template <class data_t>
    void compute(const Cell<data_t> &cell, const double y, data_t &out_re,
                 data_t &out_im) const
    {
        const double h = m_dx;
        const double yhat =
            y * axis_weight_factor(y, h, AxisWeight::flux_form);
        static const double cf[4] = {1., -27., 27., -1.}; // staggered div
        static const double it[4] = {-1. / 16., 9. / 16., 9. / 16.,
                                     -1. / 16.}; // 4-point interpolation
        static const double c1[5] = {1., -8., 0., 8., -1.}; // centred deriv

        data_t T_re = 0., T_im = 0.;

        // faces s + 1/2 (between offsets s and s + 1), s = -2 .. 1
        for (int q = 0; q < 4; ++q)
        {
            const int s = q - 2;
            const double yf = std::fabs(y + (s + 0.5) * h); // y-face position

            // x fluxes:  yhat B^xx Dx phi
            const data_t ax = yhat * face_B(cell, comp_xx, 0, 0, 0, s);
            T_re += cf[q] * ax * face_D(cell, c_phi, 0, 0, 0, s) / (24. * h);
            T_im += cf[q] * ax * face_D(cell, c_phi_Im, 0, 0, 0, s) / (24. * h);

            // y fluxes:  |y_f| (B^yy Dy phi + B^yx Iy dx phi)
            const data_t byy = yf * face_B(cell, comp_yy, 1, 0, 0, s);
            const data_t byx = yf * face_B(cell, comp_xy, 1, 0, 0, s);
            T_re += cf[q] * (byy * face_D(cell, c_phi, 1, 0, 0, s) / (24. * h) +
                             byx * face_I_dx(cell, c_phi, s));
            T_im +=
                cf[q] * (byy * face_D(cell, c_phi_Im, 1, 0, 0, s) / (24. * h) +
                         byx * face_I_dx(cell, c_phi_Im, s));
        }
        T_re /= 24. * h;
        T_im /= 24. * h;

        // transpose of the y flux of d_x phi:  dx Iy^T (|y_f| B^yx Dy phi),
        // Iy^T P at column ox = sum_s it[s] P(face s + 1/2 of column ox)
        data_t S_re = 0., S_im = 0.;
        for (int m = 0; m < 5; ++m)
        {
            const int ox = m - 2;
            if (ox == 0)
                continue;
            data_t Q_re = 0., Q_im = 0.;
            for (int q = 0; q < 4; ++q)
            {
                const int s = q - 2;
                const double yf = std::fabs(y + (s + 0.5) * h);
                const data_t P = it[q] * yf * face_B(cell, comp_xy, 1, ox, 0, s);
                Q_re += P * face_D(cell, c_phi, 1, ox, 0, s);
                Q_im += P * face_D(cell, c_phi_Im, 1, ox, 0, s);
            }
            S_re += c1[m] * Q_re;
            S_im += c1[m] * Q_im;
        }
        T_re += S_re / (24. * h * 12. * h);
        T_im += S_im / (24. * h * 12. * h);

        // L phi = -(1/W) T,   1/W = chi^{3/2} / yhat
        const data_t chi = at(cell, c_chi, 0, 0);
        const data_t inv_W = chi * sqrt(chi) / yhat;
        out_re = -inv_W * T_re;
        out_im = -inv_W * T_im;
    }

  protected:
    const double m_dx;

    enum
    {
        comp_xx,
        comp_xy,
        comp_yy
    };

    /// value of grid component ivar at offset (ox, oy) from the current cell
    template <class data_t>
    ALWAYS_INLINE data_t at(const Cell<data_t> &cell, const int ivar,
                            const int ox, const int oy) const
    {
        const auto &bp = cell.get_box_pointers();
        const int idx = cell.get_in_index() + ox * bp.m_in_stride[0] +
                        oy * bp.m_in_stride[1];
        return SIMDIFY<data_t>(bp.m_in_ptr[ivar])[idx];
    }

    /// B^ab = alpha chi^{-1/2} h^ab at offset (ox, oy)
    template <class data_t>
    ALWAYS_INLINE data_t B(const Cell<data_t> &cell, const int comp,
                           const int ox, const int oy) const
    {
        const data_t h11 = at(cell, c_h11, ox, oy);
        const data_t h12 = at(cell, c_h12, ox, oy);
        const data_t h22 = at(cell, c_h22, ox, oy);
        const data_t det = h11 * h22 - h12 * h12;
        const data_t hUU = (comp == comp_xx) ? h22 / det
                                          : ((comp == comp_xy) ? -h12 / det
                                                            : h11 / det);
        return at(cell, c_lapse, ox, oy) * hUU / sqrt(at(cell, c_chi, ox, oy));
    }

    /// offset of the o-th cell along dir from (bx, by)
    static ALWAYS_INLINE int offx(int dir, int bx, int o)
    {
        return dir == 0 ? bx + o : bx;
    }
    static ALWAYS_INLINE int offy(int dir, int by, int o)
    {
        return dir == 1 ? by + o : by;
    }

    /// 4-point interpolation of B^comp to the face between cells s and s + 1
    /// along dir, counted from the cell at offset (bx, by)
    template <class data_t>
    ALWAYS_INLINE data_t face_B(const Cell<data_t> &cell, const int comp,
                                const int dir, const int bx, const int by,
                                const int s) const
    {
        data_t r = 0.;
        static const double w[4] = {-1., 9., 9., -1.};
        for (int k = 0; k < 4; ++k)
        {
            const int o = s - 1 + k;
            r += w[k] * B(cell, comp, offx(dir, bx, o), offy(dir, by, o));
        }
        return r / 16.;
    }

    /// staggered derivative along dir at the same face, times 24 dx
    template <class data_t>
    ALWAYS_INLINE data_t face_D(const Cell<data_t> &cell, const int ivar,
                                const int dir, const int bx, const int by,
                                const int s) const
    {
        data_t r = 0.;
        static const double w[4] = {1., -27., 27., -1.};
        for (int k = 0; k < 4; ++k)
        {
            const int o = s - 1 + k;
            r += w[k] * at(cell, ivar, offx(dir, bx, o), offy(dir, by, o));
        }
        return r;
    }

    /// collocated 4th-order d_x of ivar at offset (0, oy), interpolated to
    /// the y face between oy = s and s + 1
    template <class data_t>
    ALWAYS_INLINE data_t face_I_dx(const Cell<data_t> &cell, const int ivar,
                                   const int s) const
    {
        static const double w[4] = {-1., 9., 9., -1.};
        data_t r = 0.;
        for (int k = 0; k < 4; ++k)
        {
            const int oy = s - 1 + k;
            const data_t dxv =
                (at(cell, ivar, -2, oy) - 8. * at(cell, ivar, -1, oy) +
                 8. * at(cell, ivar, 1, oy) - at(cell, ivar, 2, oy)) /
                (12. * m_dx);
            r += w[k] * dxv;
        }
        return r / 16.;
    }
};

#endif /* FLUXFORMKGCARTOON_HPP_ */
