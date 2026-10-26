// Single-box unit test of the assembled 2D Dirac code (no time evolution, runs in about a second on one core):
// 2D initial data for a star -> one evaluation of the CCZ4CartoonDirac RHS and of ConstraintsCartoonDirac, at
// dx = 0.1 and 0.05. For an unperturbed star d_t psi = -i omega psi and the metric RHS and Ham vanish up to truncation
// error, which should fall by ~16 between the two resolutions (4th order). For a perturbed star only Ham -> 0.
// Modes (third argument): 0 star at rest (default); 1 star boosted along the axis (rapidity 0.3): Ham and Mom must
// converge to zero and d_t psi must equal the exact time derivative of psi_lab = U(eta) psi_rest(t~, x~) (from
// DiracStarBoost.hpp at t = +-h); 2 equal-mass head-on binary: the data must be mirror-symmetric about x = 0, with
// psi(-x) = diag(1, -1, -1, 1) psi(x). Optional fourth argument: 1D solver gridpoints (default 20001).
// Build: make all (uses the example's UserVariables.hpp), run: ./box_test_dirac*.ex <A_central> <perturb 0|1> [mode] [n1d]
#include "BoxIterator.H"
#include "BoxLoops.hpp"
#include "CCZ4CartoonDirac.hpp"
#include "ConstraintsCartoonDirac.hpp"
#include "DiracStar.hpp"
#include "FArrayBox.H"
#include "GammaCartoonCalculator.hpp"
#include "MovingPunctureGauge.hpp"
#include "SetValue.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

int main(int argc, char **argv)
{
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif
    if (argc < 3)
    {
        printf("usage: %s <A_central> <perturb 0|1>\n", argv[0]);
        return 1;
    }
    const double A = atof(argv[1]);
    const int perturb = atoi(argv[2]);
    const int mode = (argc > 3) ? atoi(argv[3]) : 0;
    const double extent = 6.;
    DiracStar_params_t idp;
    idp.solver.A_central = A;
    idp.solver.perturb = perturb;
    idp.solver.perturb_amp = 0.01;
    idp.solver.perturb_spread = 3.;
    idp.star_centre = {0., 0.};
    idp.phase = 0.3;
    idp.write_1d_profiles = false;
    idp.rapidity = (mode == 1) ? 0.3 : 0.;
    if (argc > 4)
        idp.solver.n_gridpoints = atoi(argv[4]);
    DiracFieldCartoon::params_t dp{1.};

    CCZ4_params_t<MovingPunctureGauge::params_t> params;
    params.kappa1 = 0.1; params.kappa2 = 0.; params.kappa3 = 1.; params.covariantZ4 = true;
    params.lapse_advec_coeff = 1.; params.lapse_power = 1.; params.lapse_coeff = 2.;
    params.shift_Gamma_coeff = 0.75; params.shift_advec_coeff = 0.; params.eta = 1.;

    if (mode == 2)
    {
        idp.binary = true;
        idp.solver2 = idp.solver;
        idp.rapidity = idp.rapidity2 = 0.2;
        idp.phase = idp.phase2 = 0.25;
        idp.separation = 4.;
        idp.radius_width1 = idp.radius_width2 = 100.;
        const double dx = 0.1;
        const int nx = 60;
        Box box(IntVect(-nx, 0), IntVect(nx - 1, nx - 1));
        for (int id = 0; id < 3; ++id)
        {
            idp.id_choice = id;
            FArrayBox in(box, NUM_VARS);
            DiracStar star(idp, dx, 0);
            star.compute_1d_solution();
            BoxLoops::loop(make_compute_pack(SetValue(0.), star), in, in, box, disable_simd());
            const double par[4] = {1., -1., -1., 1.};
            double e_par = 0., e_metric = 0., maxpsi = 0.;
            for (BoxIterator bit(box); bit.ok(); ++bit)
            {
                IntVect iv = bit(), ivm = iv;
                ivm[0] = -iv[0] - 1;
                for (int c = 0; c < 4; ++c)
                    for (int ri = 0; ri < 2; ++ri)
                    {
                        const int k = c_psi0_Re + 2 * c + ri;
                        e_par = std::max(e_par, std::abs(in(ivm, k) - par[c] * in(iv, k)));
                        maxpsi = std::max(maxpsi, std::abs(in(iv, k)));
                    }
                e_metric = std::max(e_metric, std::abs(in(ivm, c_chi) - in(iv, c_chi)) +
                                                  std::abs(in(ivm, c_lapse) - in(iv, c_lapse)) +
                                                  std::abs(in(ivm, c_shift1) + in(iv, c_shift1)) +
                                                  std::abs(in(ivm, c_K) - in(iv, c_K)) +
                                                  std::abs(in(ivm, c_h12) + in(iv, c_h12)));
            }
            printf("binary id_choice %d: max |psi(-x) - P psi(x)| = %.1e (max |psi| %.2e), metric/gauge mirror error "
                   "%.1e\n", id, e_par, maxpsi, e_metric);
        }
#ifdef CH_MPI
        MPI_Finalize();
#endif
        return 0;
    }

    for (double dx : {0.1, 0.05})
    {
        const int nx = static_cast<int>(round(extent / dx));
        Box box(IntVect(-nx, 0), IntVect(nx - 1, nx - 1));
        Box gbox = box; gbox.grow(6);
        Box mbox = box; mbox.grow(3);
        FArrayBox in(gbox, NUM_VARS), rhs(box, NUM_VARS), diag(box, NUM_DIAGNOSTIC_VARS);
        DiracStar star(idp, dx, 0);
        star.compute_1d_solution();
        BoxLoops::loop(make_compute_pack(SetValue(0.), star), in, in, gbox, disable_simd());
        BoxLoops::loop(GammaCartoonCalculator(dx), in, in, mbox, disable_simd());
        BoxLoops::loop(CCZ4CartoonDirac<MovingPunctureGauge, FourthOrderDerivatives>(params, dx, 0.3, dp, 1., 0), in,
                       rhs, box);
        BoxLoops::loop(ConstraintsCartoonDirac(dx, dp, 1.), in, diag, box);

        const double w = star.m_1d_sol.omega;
        FArrayBox exact_dt(box, 8);
        if (mode == 1)
        {
            const double h = 1e-4;
            for (BoxIterator bit(box); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                const double x = (iv[0] + 0.5) * dx, y = (iv[1] + 0.5) * dx;
                double pR[4], pI[4], mR[4], mI[4];
                DiracStar::spinor(boosted_dirac_star(star.m_1d_sol, x, y, 0., idp.rapidity, idp.phase, h), pR, pI);
                DiracStar::spinor(boosted_dirac_star(star.m_1d_sol, x, y, 0., idp.rapidity, idp.phase, -h), mR, mI);
                for (int c = 0; c < 4; ++c)
                {
                    exact_dt(iv, 2 * c) = (pR[c] - mR[c]) / (2. * h);
                    exact_dt(iv, 2 * c + 1) = (pI[c] - mI[c]) / (2. * h);
                }
            }
        }
        double e_psi = 0., n_psi = 0., e_metric = 0., ham = 0., rho16 = 0., mom = 0., modpsi_c = 0., Nsum = 0.;
        for (BoxIterator bit(box); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            for (int c = 0; c < 4; ++c)
            {
                // at rest d_t psi = -i omega psi:  Re' = omega Im,  Im' = -omega Re
                const double re = in(iv, c_psi0_Re + 2 * c), im = in(iv, c_psi0_Im + 2 * c);
                const double dre = (mode == 1) ? exact_dt(iv, 2 * c) : w * im;
                const double dim = (mode == 1) ? exact_dt(iv, 2 * c + 1) : -w * re;
                e_psi = std::max(e_psi, std::max(std::abs(rhs(iv, c_psi0_Re + 2 * c) - dre),
                                                 std::abs(rhs(iv, c_psi0_Im + 2 * c) - dim)));
                n_psi = std::max(n_psi, std::sqrt(dre * dre + dim * dim));
            }
            for (int c = 0; c < c_psi0_Re; ++c)
                e_metric = std::max(e_metric, std::abs(rhs(iv, c)));
            ham = std::max(ham, std::abs(diag(iv, c_Ham)));
            rho16 = std::max(rho16, 16. * M_PI * std::abs(diag(iv, c_rho)));
            mom = std::max(mom, std::abs(diag(iv, c_Mom)));
            modpsi_c = std::max(modpsi_c, diag(iv, c_mod_psi));
            Nsum += diag(iv, c_N) * dx * dx;
        }
        if (mode == 1)
            e_metric = 0.; // a boosted star moves
        printf("mode %d A = %.2f perturb %d dx = %.3f: |d_t psi - exact|/|exact| = %.2e  max|metric rhs| = %.2e  "
               "max|Ham|/max(16 pi rho) = %.2e  max|Mom| = %.2e  max mod_psi = %.6f (A_central %.6f)  N(|x|,y<6) = %.5f\n",
               mode, A, perturb, dx, e_psi / n_psi, e_metric, ham / rho16, mom, modpsi_c, star.m_1d_sol.A_central, Nsum);
    }
#ifdef CH_MPI
    MPI_Finalize();
#endif
    return 0;
}
