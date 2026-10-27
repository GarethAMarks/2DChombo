// Single-box unit test of the assembled 2D Proca code (no time evolution, runs in about a second on one core):
// 2D initial data for a star -> one evaluation of the CCZ4CartoonProca RHS and of ConstraintsCartoonProca, at
// dx = 0.1 and 0.05 on [-6, 6] x [0, 6] (the axis included; the ghost cells at y < 0 come from the initial data). For an
// unperturbed star d_t (X_i, E^i, X_phi) = -i omega (X_i, E^i, X_phi), d_t Z = 0, and the metric RHS, Ham and the Gauss
// constraint vanish up to truncation error, which should fall by ~16 between the two resolutions (4th order). For a
// perturbed star only the constraints -> 0.
// Modes (third argument): 0 star at rest (default); 1 star boosted along the axis (rapidity 0.3): the constraints must
// converge to zero and the Proca RHS must equal the exact time derivative of the boosted solution (from
// ProcaStarBoost.hpp at t = +-h); 2 equal-mass head-on binary: the data must be mirror-symmetric about x = 0, with
// X_x and E^x odd; 3 aspherical (axisymmetric) perturbation m = 3, delta = 0.05 of a star at rest (arXiv:2609.30021
// Eq. (11)): the Proca variables must be F = exp(-2 delta cos(m phi)) times the unperturbed ones (phi the polar angle
// from the axis), the Gauss constraint and Mom must still converge to zero (E^i is radial and F depends only on angles),
// while Ham is O(delta) (the constraints are not re-solved). Optional fourth argument: 1D solver gridpoints (default
// 30001).
// Build: make all (uses the example's UserVariables.hpp), run: ./box_test_proca*.ex <A_central> <perturb 0|1> [mode] [n1d]
#include "BoxIterator.H"
#include "BoxLoops.hpp"
#include "CCZ4CartoonProca.hpp"
#include "ConstraintsCartoonProca.hpp"
#include "FArrayBox.H"
#include "GammaCartoonCalculator.hpp"
#include "MovingPunctureGauge.hpp"
#include "ProcaStar.hpp"
#include "SetValue.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv)
{
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif
    if (argc < 3)
    {
        printf("usage: %s <A_central> <perturb 0|1> [mode] [n1d]\n", argv[0]);
        return 1;
    }
    const double A = atof(argv[1]);
    const int perturb = atoi(argv[2]);
    const int mode = (argc > 3) ? atoi(argv[3]) : 0;
    const double extent = 6.;
    ProcaStar_params_t idp;
    idp.solver.A_central = A;
    idp.solver.perturb = perturb;
    idp.solver.perturb_amp = 0.01;
    idp.solver.perturb_spread = 3.;
    idp.star_centre = {0., 0.};
    idp.phase = 0.3;
    idp.write_1d_profiles = false;
    idp.rapidity = (mode == 1) ? 0.3 : 0.;
    if (mode == 3)
    {
        idp.aspherical.active = true;
        idp.aspherical.m = 3;
        idp.aspherical.delta = 0.05;
    }
    if (argc > 4)
        idp.solver.n_gridpoints = atoi(argv[4]);
    ProcaFieldCartoon::params_t pp{1., 1.};

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
            ProcaStar star(idp, dx, 0);
            star.compute_1d_solution();
            BoxLoops::loop(make_compute_pack(SetValue(0.), star), in, in, box, disable_simd());
            double e_par = 0., e_metric = 0., maxP = 0.;
            for (BoxIterator bit(box); bit.ok(); ++bit)
            {
                IntVect iv = bit(), ivm = iv;
                ivm[0] = -iv[0] - 1;
                for (int k = c_Xx_Re; k < NUM_VARS; ++k)
                {
                    const bool odd = (k == c_Xx_Re || k == c_Xx_Im || k == c_Ex_Re || k == c_Ex_Im);
                    e_par = std::max(e_par, std::abs(in(ivm, k) - (odd ? -1. : 1.) * in(iv, k)));
                    maxP = std::max(maxP, std::abs(in(iv, k)));
                }
                e_metric = std::max(e_metric, std::abs(in(ivm, c_chi) - in(iv, c_chi)) +
                                                  std::abs(in(ivm, c_lapse) - in(iv, c_lapse)) +
                                                  std::abs(in(ivm, c_shift1) + in(iv, c_shift1)) +
                                                  std::abs(in(ivm, c_K) - in(iv, c_K)) +
                                                  std::abs(in(ivm, c_h12) + in(iv, c_h12)));
            }
            printf("binary id_choice %d: max |Proca var(-x) - parity var(x)| = %.1e (max |var| %.2e), metric/gauge "
                   "mirror error %.1e\n", id, e_par, maxP, e_metric);
        }
#ifdef CH_MPI
        MPI_Finalize();
#endif
        return 0;
    }

    const int nP = NUM_VARS - c_Xx_Re; // 12 Proca variables
    for (double dx : {0.1, 0.05})
    {
        const int nx = static_cast<int>(round(extent / dx));
        Box box(IntVect(-nx, 0), IntVect(nx - 1, nx - 1));
        Box gbox = box; gbox.grow(6);
        Box mbox = box; mbox.grow(3);
        FArrayBox in(gbox, NUM_VARS), rhs(box, NUM_VARS), diag(box, NUM_DIAGNOSTIC_VARS);
        ProcaStar star(idp, dx, 0);
        star.compute_1d_solution();
        const ProcaStarSolution &sol = star.m_1d_sol;
        BoxLoops::loop(make_compute_pack(SetValue(0.), star), in, in, gbox, disable_simd());
        BoxLoops::loop(GammaCartoonCalculator(dx), in, in, mbox, disable_simd());
        double e_factor = 0., F_min = 1., F_max = 1.;
        if (mode == 3)
        {
            // the perturbed Proca variables must be F times the unperturbed ones, the rest unchanged
            ProcaStar_params_t ref = idp;
            ref.aspherical.active = false;
            ProcaStar star_ref(ref, dx, 0);
            star_ref.compute_1d_solution();
            FArrayBox in_ref(box, NUM_VARS);
            BoxLoops::loop(make_compute_pack(SetValue(0.), star_ref), in_ref, in_ref, box, disable_simd());
            for (BoxIterator bit(box); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                const double F = idp.aspherical.factor((iv[0] + 0.5) * dx, (iv[1] + 0.5) * dx, 0.);
                F_min = std::min(F_min, F);
                F_max = std::max(F_max, F);
                for (int c = 0; c < NUM_VARS; ++c)
                    e_factor = std::max(e_factor, std::abs(in(iv, c) - ((c >= c_Xx_Re) ? F : 1.) * in_ref(iv, c)));
            }
        }
        BoxLoops::loop(CCZ4CartoonProca<MovingPunctureGauge, FourthOrderDerivatives>(params, dx, 0.3, pp, 1., 0), in,
                       rhs, box);
        BoxLoops::loop(ConstraintsCartoonProca(dx, pp, 1.), in, diag, box);

        const double w = sol.omega;
        double gauss_out = 0., e_P = 0., n_P = 0., e_metric = 0., ham = 0., rho16 = 0., mom = 0., gauss = 0., mu2X = 0., e_n = 0.,
               n_max = 0., modX = 0.;
        for (BoxIterator bit(box); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            const double x = (iv[0] + 0.5) * dx, y = (iv[1] + 0.5) * dx;
            double exact[12];
            if (mode == 1)
            {
                const double h = 1e-4;
                const auto p = boosted_proca_star(sol, x, y, 0., idp.rapidity, idp.phase, h);
                const auto m = boosted_proca_star(sol, x, y, 0., idp.rapidity, idp.phase, -h);
                auto d = [&](std::complex<double> a, std::complex<double> b) { return (a - b) / (2. * h); };
                for (int i = 0; i < 2; ++i)
                {
                    exact[i] = d(p.X[i], m.X[i]).real();
                    exact[2 + i] = d(p.X[i], m.X[i]).imag();
                    exact[4 + i] = d(p.E[i], m.E[i]).real();
                    exact[6 + i] = d(p.E[i], m.E[i]).imag();
                }
                exact[8] = d(p.Xphi, m.Xphi).real();
                exact[9] = d(p.Xphi, m.Xphi).imag();
            }
            else
            {
                // -i omega: d_t Re = omega Im, d_t Im = -omega Re
                for (int i = 0; i < 2; ++i)
                {
                    exact[i] = w * in(iv, c_Xx_Im + i);
                    exact[2 + i] = -w * in(iv, c_Xx_Re + i);
                    exact[4 + i] = w * in(iv, c_Ex_Im + i);
                    exact[6 + i] = -w * in(iv, c_Ex_Re + i);
                }
                exact[8] = w * in(iv, c_Xphi_Im);
                exact[9] = -w * in(iv, c_Xphi_Re);
            }
            exact[10] = exact[11] = 0.;
            for (int k = 0; k < nP; ++k)
            {
                e_P = std::max(e_P, std::abs(rhs(iv, c_Xx_Re + k) - exact[k]));
                n_P = std::max(n_P, std::abs(exact[k]));
            }
            for (int c = 0; c < c_Xx_Re; ++c)
                e_metric = std::max(e_metric, std::abs(rhs(iv, c)));
            ham = std::max(ham, std::abs(diag(iv, c_Ham)));
            rho16 = std::max(rho16, 16. * M_PI * std::abs(diag(iv, c_rho)));
            mom = std::max(mom, std::abs(diag(iv, c_Mom)));
            gauss = std::max(gauss, diag(iv, c_Gauss));
            if (std::hypot(x, y) > 0.5)
                gauss_out = std::max(gauss_out, diag(iv, c_Gauss)); // away from the centre, where F is undefined
            mu2X = std::max(mu2X, std::hypot(in(iv, c_Xphi_Re), in(iv, c_Xphi_Im)));
            modX = std::max(modX, diag(iv, c_mod_Xphi));
            if (mode == 0)
            {
                const double r = std::hypot(x, y), psi = sol.get_psi(r);
                const double n_exact = -sol.get_g(r) * sol.get_E(r) / (sol.get_lapse(r) * pow(psi, 4));
                e_n = std::max(e_n, std::abs(diag(iv, c_n) - n_exact));
                n_max = std::max(n_max, std::abs(n_exact));
            }
        }
        if (mode == 1)
            e_metric = 0.; // a boosted star moves
        if (mode == 3)
        {
            printf("mode 3 (aspherical m = 3, delta = 0.05) dx = %.3f: max |vars - F vars_unperturbed| = %.1e "
                   "(F in [%.4f, %.4f])  max|Gauss|/max|mu^2 X_phi| = %.2e (%.2e for r > 0.5)  max|Mom| = %.2e  "
                   "max|Ham|/max(16 pi rho) = %.2e\n",
                   dx, e_factor, F_min, F_max, gauss / mu2X, gauss_out / mu2X, mom, ham / rho16);
            continue;
        }
        printf("mode %d A = %.3f perturb %d dx = %.3f: |d_t Proca - exact|/|exact| = %.2e  max|metric rhs| = %.2e  "
               "max|Ham|/max(16 pi rho) = %.2e  max|Mom| = %.2e  max|Gauss|/max|mu^2 X_phi| = %.2e\n",
               mode, A, perturb, dx, e_P / n_P, e_metric, ham / rho16, mom, gauss / mu2X);
        printf("   max mod_Xphi = %.6f (A_central %.6f)", modX, sol.A_central);
        if (mode == 0)
            printf("   max |n - n_exact| / max n_exact = %.1e", e_n / n_max);
        printf("\n");
    }
#ifdef CH_MPI
    MPI_Finalize();
#endif
    return 0;
}
