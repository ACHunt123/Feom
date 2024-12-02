// File: general.h
#ifndef GENERAL_H
#define GENERAL_H

#include <mpi.h>
#include <complex>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#define ARMA_ALLOW_FAKE_GCC
//#define ARMA_NO_DEBUG
#include <armadillo>
namespace a = arma;

using Int = long;
using SuperInt = unsigned long long int;
using Real = double;
using Cplx = std::complex<double>;
using String = std::string;

using Vec = a::Col<Real>;
using IntVec = a::Col<Int>;
using CplxVec = a::Col<Cplx>;
using Mat = a::Mat<Real>;
using IntMat = a::Mat<Int>;
using CplxMat = a::Mat<Cplx>;
using Cube = a::Cube<Real>;
using CplxCube = a::Cube<Cplx>;
template<class T>
using VecGen = a::Col<T>; 
template<class T>
using CppVec = std::vector<T>;

extern const Real pi;
extern const Real hbar;
extern const Real au2fs;
extern const Real fs2au;
extern const Cplx imu;

enum Heom_switch
{
    HEOM_NONE,
    HEOM_POTENTIAL,
    HEOM_TCF_QQ,
    HEOM_KUBO_QQ,
    HEOM_TCF_Q2Q2,
    HEOM_KUBO_Q2Q2,
};

const std::map<String,Heom_switch> heom_switch_map =
{
    {"none", HEOM_NONE},
    {"potential", HEOM_POTENTIAL,},
    {"tcf_qq", HEOM_TCF_QQ,   },
    {"kubo_qq", HEOM_KUBO_QQ,  },
    {"tcf_q2q2", HEOM_TCF_Q2Q2, },
    {"kubo_q2q2", HEOM_KUBO_Q2Q2,},
};

enum Bath_type
{
    NO_BATH,
    WHITE,
    DEBYE,
    DEBYE_CORRECTION,
    DEBYE_CORRECTION2,
};
const std::map<String,Bath_type> bath_type_map =
{
    {"no_bath", NO_BATH},
    {"white", WHITE},
    {"debye", DEBYE}, // Do not use
    {"debye_correction", DEBYE_CORRECTION},
    {"debye_correction2", DEBYE_CORRECTION2},
};

// File naming functions
String get_file_real_au_dat(String name, String kind);
String get_file_real_au_err(String name, String kind);
String get_file_imag_au_dat(String name, String kind);
String get_file_imag_au_err(String name, String kind);
String get_file_real_fs_dat(String name, String kind);
String get_file_real_fs_err(String name, String kind);
String get_file_imag_fs_dat(String name, String kind);
String get_file_imag_fs_err(String name, String kind);

SuperInt factorial(Int N);
SuperInt binom(Int n, Int k);
IntMat get_flattenting_coefficients(Int K, Int L);

#endif
