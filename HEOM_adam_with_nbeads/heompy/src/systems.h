#ifndef SYSTEMS_H
#define SYSTEMS_H
#include "general.h"
#include "input.h"
class SysBase
{
public:
    Real t0{},t1{};
    Real dt{};
    Real dt_sample{};
    Int n_t{},tau{},t_offset{};
    Vec ts;

    Int sample_step{};
    Bath_type bath_type;
    Heom_switch heom_switch;
    Real eta{};
    Real omega_c{};
    Real beta{};
    Real mass{};
    Int cL{};
    Int cK{};
    Int n_ee{};

    CplxVec cs;
    Vec gammas;
    Mat q_mat, dq_mat, qq_mat;
    CplxMat rho0, ham, white_ham, tmp_ham;
    IntMat flattening;
    IntMat indices;
    Int n_ados{};
    CplxCube ados, k1, k2, k3, k4;

    int rank{};
    int n_procs{};
    IntMat mpi_split;

    Cplx tmpcplx{};
    bool stability{};
    CplxMat rho, iham, q_rho, rho_q, comm, iq_mat;
    Cplx low_T_coef_R{}, low_T_coef_I{};
    IntVec curvec{}, nextvec{};
    Int curind{}, nextind{};

    SysBase(const Input &inp);
    Int z(const IntVec & intvec);
    Real get_magnitude();
    void write_tcf_qq(const String & filename, const Real & time);
    void write_tcf_q2q2(const String & filename, const Real & time);
    void write_tcf_q4q4(const String & filename, const Real & time);
    bool check_stability(const Real & maxval);
    void postmultiply_by_qs();
    void derivative(const CplxCube & ado, const Real & coef, CplxCube & drdt);
    void propagate();
    void change_to_white_ham();
    void change_from_white_ham();
};
#endif
