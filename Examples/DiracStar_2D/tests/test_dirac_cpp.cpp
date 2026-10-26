// Standalone test of Source/Matter/DiracFieldCartoon.hpp (no Chombo): reads one point per line from stdin
//   chi h00 h01 h11 hww K A00 A01 A11 Aww lapse b0 b1
//   [d_x, d_y of] chi h00 h01 h11 hww lapse b0 b1
//   psiRe[4] psiIm[4]  dpsiRe[4][x,y] dpsiIm[4][x,y]  y mu
// and prints the 8 RHS components (no advection) and rho Sx Sy Sxx Sxy Syy Sww n.
#include <cstdio>
#include <iostream>
enum
{
    c_psi0_Re, c_psi0_Im, c_psi1_Re, c_psi1_Im, c_psi2_Re, c_psi2_Im, c_psi3_Re, c_psi3_Im
};
#include "DiracFieldCartoon.hpp"

template <class data_t> struct TestVars : DiracSpinorVars<data_t>
{
    data_t chi, hww, K, Aww, lapse;
    Tensor<2, data_t> h, A;
    Tensor<1, data_t> shift;
};

int main()
{
    TestVars<double> v;
    TestVars<Tensor<1, double>> d1;
    double y, mu;
    while (std::cin >> v.chi >> v.h[0][0] >> v.h[0][1] >> v.h[1][1] >> v.hww >> v.K >> v.A[0][0] >> v.A[0][1] >>
           v.A[1][1] >> v.Aww >> v.lapse >> v.shift[0] >> v.shift[1])
    {
        v.h[1][0] = v.h[0][1];
        v.A[1][0] = v.A[0][1];
        for (int k = 0; k < 2; ++k)
        {
            std::cin >> d1.chi[k] >> d1.h[0][0][k] >> d1.h[0][1][k] >> d1.h[1][1][k] >> d1.hww[k] >> d1.lapse[k] >>
                d1.shift[0][k] >> d1.shift[1][k];
            d1.h[1][0][k] = d1.h[0][1][k];
        }
        for (int c = 0; c < 4; ++c)
            std::cin >> v.psi_Re[c];
        for (int c = 0; c < 4; ++c)
            std::cin >> v.psi_Im[c];
        for (int c = 0; c < 4; ++c)
            for (int k = 0; k < 2; ++k)
                std::cin >> d1.psi_Re[c][k];
        for (int c = 0; c < 4; ++c)
            for (int k = 0; k < 2; ++k)
                std::cin >> d1.psi_Im[c][k];
        std::cin >> y >> mu;
        DiracFieldCartoon dirac({mu});
        DiracSpinorVars<double> rhs;
        DiracSources_t<double> src;
        dirac.compute(rhs, src, v, d1, y);
        for (int c = 0; c < 4; ++c)
            printf("%.17e %.17e ", rhs.psi_Re[c], rhs.psi_Im[c]);
        printf("%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", src.rho, src.Si[0], src.Si[1], src.Sij[0][0],
               src.Sij[0][1], src.Sij[1][1], src.Sww, src.n);
    }
    return 0;
}
