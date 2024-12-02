// File: general.cpp
#include "general.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib> // for exit()

// Constants
extern const Real pi{3.141592653589793238462643383279502884197169399375105L};
extern const Real hbar{1.000000000000L};
extern const Real au2fs{2.418884326509e-2L};
extern const Real fs2au{1/2.418884326509e-2L};
extern const Cplx imu{0.0,1.0};

// File naming functions

String get_file_real_au_dat(String name, String kind)
{
    return name + "_" + kind + "_real_au_dat";
}

String get_file_real_au_err(String name, String kind)
{
    return name + "_" + kind + "_real_au_err";
}

String get_file_imag_au_dat(String name, String kind)
{
    return name + "_" + kind + "_imag_au_dat";
}

String get_file_imag_au_err(String name, String kind)
{
    return name + "_" + kind + "_imag_au_err";
}

String get_file_real_fs_dat(String name, String kind)
{
    return name + "_" + kind + "_real_fs_dat";
}

String get_file_real_fs_err(String name, String kind)
{
    return name + "_" + kind + "_real_fs_err";
}

String get_file_imag_fs_dat(String name, String kind)
{
    return name + "_" + kind + "_imag_fs_dat";
}

String get_file_imag_fs_err(String name, String kind)
{
    return name + "_" + kind + "_imag_fs_err";
}

// Factorial

SuperInt factorial(Int N)
{
    SuperInt result{1};
    if (N>0)
    {
        for (Int i=N; i>0; i=i-1)
        {
            result = result*i;
        }
    }
    else if (N==0)
    {
        result = 1;
    }
    else
    {
        std::cerr << "Factorial is not defined for negative numbers.\n";
        std::exit(1);
    }
    return result;
}
//SuperInt binom(Int n, Int k)
//{
//    SuperInt result{1};
//    if (n>=k && n>=0 && k>=0)
//    {
//        Int lower{};
//        if ((n-k)<k) lower = n-k;
//        else lower = k;
//        for(Int i=n; i>(n-lower); i=i-1)
//        {
//            result = result*i;
//        }
//        return result/factorial(lower);
//    }
//    else
//    {
//        std::cerr << "Invalid binomial input.\n";
//        std::exit(1);
//
//    }
//}

SuperInt binom(Int n, Int k)
{
   if (k == 0 || k == n)
   return 1;
   return binom(n - 1, k - 1) + binom(n - 1, k);
}

IntMat get_flattenting_coefficients(Int K, Int L)
{
    IntMat mat(L+1,K+1,a::fill::zeros);
    for (int k=0; k<=K; k++)
    {
        mat(0,k) = 0;
        for (int n=1; n<=L; n++)
        {
            mat(n,k) = binom(n+k,k+1);
        }
    }
    return mat;
}
