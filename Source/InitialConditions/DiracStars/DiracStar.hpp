/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIRACSTAR_HPP_
#define DIRACSTAR_HPP_

#include "CCZ4CartoonDirac.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DiracStarSolution.hpp"
#include "MayDay.H"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <array>
#include <cmath>

//! Parameters of a single Dirac star
struct DiracStar_params_t
{
    DiracStarSolution::params_t solver;      //!< 1D solver parameters (mu, A_central, eigen, perturbation, ...)
    std::array<double, CH_SPACEDIM> star_centre; //!< centre of the star; must lie on the axis (y = 0)
    double phase;                            //!< overall phase of the spinor, in units of pi
    bool write_1d_profiles;                  //!< write the 1D solution to DiracStar_1d_profiles.dat (rank 0)
};

/// Initial data for a Dirac star in the 2D cartoon code
/**
 * Solves the spherical star with DiracStarSolution (a port of the Spherical_BS_Evolver solver, optionally perturbed)
 * and sets, in isotropic coordinates centred on star_centre (on the symmetry axis),
 *   chi = psi^{-4}, h_ij = delta_ij, h_ww = 1, lapse = alpha, K = A_ij = Theta = shift = B = 0,
 *   psi0 = g e^{i theta} / sqrt(4 pi),  psi1 = 0,
 *   psi2 = -i f (x/r) e^{i theta} / sqrt(4 pi),  psi3 = -i f (y/r) e^{i theta} / sqrt(4 pi),
 * the J_x = +1/2 spinor of the star (the J_x = -1/2 one is its mirror image, see DiracFieldCartoon.hpp), theta = phase pi.
 * Gamma^i must be computed afterwards (GammaCartoonCalculator). Use with disable_simd() as for BosonStar.
 **/
class DiracStar
{
  public:
    DiracStar(DiracStar_params_t a_params, double a_dx, int a_verbosity)
        : m_params(a_params), m_dx(a_dx), m_verbosity(a_verbosity)
    {
    }

    //! Computes the 1d solution (every rank solves it; it is cheap)
    void compute_1d_solution()
    {
        if (std::abs(m_params.star_centre[1]) > 1e-12)
            MayDay::Error("DiracStar: star_centre must lie on the symmetry axis (y = 0) in the cartoon code");
        if (!m_1d_sol.solve(m_params.solver, pout(), m_verbosity == 0))
            MayDay::Error("DiracStar: the 1D Dirac star solver failed (see pout); try a different "
                          "dirac_frequency_guess / dirac_alpha_central_guess or a larger dirac_solver_radius");
        pout() << "DiracStar: A_central = " << m_1d_sol.A_central << ", omega = " << m_1d_sol.omega
               << ", M = " << m_1d_sol.M << ", N = " << m_1d_sol.noether_charge
               << ", binding energy = " << m_1d_sol.binding_energy << ", r_99 = " << m_1d_sol.r_99 << endl;
#ifdef CH_MPI
        int rank;
        MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
        int rank = 0;
#endif
        if (m_params.write_1d_profiles && rank == 0)
            m_1d_sol.write_profiles("DiracStar_1d_profiles.dat");
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        CCZ4CartoonDiracVars::VarsWithGauge<data_t> vars;
        VarsTools::assign(vars, 0.);

        Coordinates<data_t> coords(current_cell, m_dx, m_params.star_centre);
        const data_t x = coords.x;
        const double y = coords.y;
        const data_t r = sqrt(x * x + y * y);

        const data_t g = m_1d_sol.get_g(r);
        const data_t f = m_1d_sol.get_f(r);
        const data_t psi_conf = m_1d_sol.get_psi(r);

        vars.chi = pow(psi_conf, -4.);
        vars.h[0][0] = 1.;
        vars.h[1][1] = 1.;
        vars.hww = 1.;
        vars.lapse = m_1d_sol.get_lapse(r);

        const double norm = 1. / sqrt(4. * M_PI);
        const double th = m_params.phase * M_PI;
        const double c = cos(th), s = sin(th);
        // psi0 = g e^{i th}; psi2,3 = -i f n e^{i th} = f n (sin th - i cos th)
        vars.psi_Re[0] = norm * g * c;
        vars.psi_Im[0] = norm * g * s;
        vars.psi_Re[2] = norm * f * (x / r) * s;
        vars.psi_Im[2] = -norm * f * (x / r) * c;
        vars.psi_Re[3] = norm * f * (y / r) * s;
        vars.psi_Im[3] = -norm * f * (y / r) * c;

        current_cell.store_vars(vars);
    }

    DiracStarSolution m_1d_sol;

  protected:
    DiracStar_params_t m_params;
    double m_dx;
    int m_verbosity;
};

#endif /* DIRACSTAR_HPP_ */
