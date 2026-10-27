// Standalone test of Source/Matter/ProcaFieldCartoon.hpp (no Chombo): reads one point per line from stdin
//   chi h00 h01 h11 hww K lapse b0 b1
//   [d_x, d_y of] chi h00 h01 h11 hww lapse b0 b1
//   X_Re[2] X_Im[2] E_Re[2] E_Im[2] Xphi_Re Xphi_Im Z_Re Z_Im,  [d_x, d_y of] these 12,
//   [d_j d_k for (j, k) = xx, xy, yx, yy of] X_Re[2] X_Im[2],   y mu kappa
// and prints the 12 RHS components (no advection, in the input order), rho Sx Sy Sxx Sxy Syy Sww n C_Re C_Im.
#include <cstdio>
#include <iostream>
enum
{
    c_Xx_Re, c_Xy_Re, c_Xx_Im, c_Xy_Im, c_Ex_Re, c_Ey_Re, c_Ex_Im, c_Ey_Im, c_Xphi_Re, c_Xphi_Im, c_Z_Re, c_Z_Im
};
#include "ProcaFieldCartoon.hpp"

template <class data_t> struct TestVars : ProcaCartoonVars<data_t>
{
    data_t chi, hww, K, lapse;
    Tensor<2, data_t> h;
    Tensor<1, data_t> shift;
};

template <class T> void read12(T *X_Re, T *X_Im, T *E_Re, T *E_Im, T &a, T &b, T &c, T &d)
{
    std::cin >> X_Re[0] >> X_Re[1] >> X_Im[0] >> X_Im[1] >> E_Re[0] >> E_Re[1] >> E_Im[0] >> E_Im[1] >> a >> b >> c >> d;
}

int main()
{
    TestVars<double> v;
    TestVars<Tensor<1, double>> d1;
    ProcaCartoonDiff2Vars<Tensor<2, double>> d2;
    double y, mu, kappa;
    while (std::cin >> v.chi >> v.h[0][0] >> v.h[0][1] >> v.h[1][1] >> v.hww >> v.K >> v.lapse >> v.shift[0] >>
           v.shift[1])
    {
        v.h[1][0] = v.h[0][1];
        for (int k = 0; k < 2; ++k)
        {
            std::cin >> d1.chi[k] >> d1.h[0][0][k] >> d1.h[0][1][k] >> d1.h[1][1][k] >> d1.hww[k] >> d1.lapse[k] >>
                d1.shift[0][k] >> d1.shift[1][k];
            d1.h[1][0][k] = d1.h[0][1][k];
        }
        read12(v.X_Re, v.X_Im, v.E_Re, v.E_Im, v.Xphi_Re, v.Xphi_Im, v.Z_Re, v.Z_Im);
        for (int k = 0; k < 2; ++k)
        {
            double s[12];
            for (double &x : s)
                std::cin >> x;
            for (int i = 0; i < 2; ++i)
            {
                d1.X_Re[i][k] = s[i];
                d1.X_Im[i][k] = s[2 + i];
                d1.E_Re[i][k] = s[4 + i];
                d1.E_Im[i][k] = s[6 + i];
            }
            d1.Xphi_Re[k] = s[8];
            d1.Xphi_Im[k] = s[9];
            d1.Z_Re[k] = s[10];
            d1.Z_Im[k] = s[11];
        }
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
            {
                for (int i = 0; i < 2; ++i)
                    std::cin >> d2.X_Re[i][j][k];
                for (int i = 0; i < 2; ++i)
                    std::cin >> d2.X_Im[i][j][k];
            }
        std::cin >> y >> mu >> kappa;
        ProcaFieldCartoon proca({mu, kappa});
        ProcaCartoonVars<double> rhs;
        ProcaSources_t<double> src;
        proca.compute(rhs, src, v, d1, d2, y);
        printf("%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e ", rhs.X_Re[0], rhs.X_Re[1],
               rhs.X_Im[0], rhs.X_Im[1], rhs.E_Re[0], rhs.E_Re[1], rhs.E_Im[0], rhs.E_Im[1], rhs.Xphi_Re, rhs.Xphi_Im,
               rhs.Z_Re, rhs.Z_Im);
        printf("%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", src.rho, src.Si[0], src.Si[1],
               src.Sij[0][0], src.Sij[0][1], src.Sij[1][1], src.Sww, src.n, src.C_Re, src.C_Im);
    }
    return 0;
}
