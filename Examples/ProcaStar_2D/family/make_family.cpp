// Regenerates proca_family.dat (ground-state mini Proca stars, mu = G = 1) with the ported solver.
//   g++ -O2 -I../../../Source/InitialConditions/ProcaStars make_family.cpp -o make_family && ./make_family > proca_family.dat
// columns: A_central M N E_b omega r_99 f(0)
#include "ProcaStarSolution.hpp"
#include <cstdio>
#include <vector>
int main()
{
    std::vector<double> As;
    for (double A = 0.001; A < 0.01; A += 0.0005) As.push_back(A);
    for (double A = 0.01; A < 0.1; A += 0.001) As.push_back(A);
    for (double A = 0.1; A < 0.6; A += 0.0025) As.push_back(A);
    for (double A = 0.61; A < 4.0; A += 0.01) As.push_back(A);
    printf("# Ground-state mini Proca stars (mu = G = 1) from ProcaStarSolution (the port of Spherical_BS_Evolver's ProcaStar)\n");
    printf("# A_central = f(0)/alpha(0)   M   N (Noether charge)   E_b = M - mu N   omega   r_99   f(0)\n");
    for (double A : As)
    {
        ProcaStarSolution::params_t p;
        p.A_central = A;
        p.alpha_central = 0.9;
        p.frequency_guess = 1.8;
        // the radius shrinks with A (r_99 = 43 at A = 0.001, 15 at 0.02, 6 at 0.6)
        p.R = (A < 0.004) ? 300. : (A < 0.02 ? 200. : (A < 0.6 ? 100. : 60.));
        p.n_gridpoints = (A < 0.004) ? 60001 : (A < 0.02 ? 40001 : (A < 0.6 ? 20001 : 24001));
        ProcaStarSolution sol;
        if (!sol.solve(p, std::cerr, true))
        {
            fprintf(stderr, "fail A=%g\n", A);
            continue;
        }
        printf("%.6f %.10f %.10f %.10e %.10f %.6f %.8e\n", A, sol.M, sol.noether_charge, sol.binding_energy, sol.omega,
               sol.r_99, sol.f0);
        fflush(stdout);
    }
}
