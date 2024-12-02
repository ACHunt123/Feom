// File: systems.cpp
#include <iostream>
#include <cstdlib>
#include <complex>
#include "general.h"
#include "utils.h"
#include "systems.h"
#include "input.h"
//****************************************************************************
// SysBase
//****************************************************************************
IntMat return_larger_indices(const IntMat & orig, const Int & L)
{
    IntMat output(orig.n_rows+1,0,a::fill::zeros);
    if (orig.n_elem>0)
    {
        for (Int i=0; i<orig.n_cols; i++)
        {
            for (Int nk=0; nk+a::sum(orig.col(i))<=L; nk++)
            {
                IntVec tmp(orig.n_rows+1);
                tmp.head(orig.n_rows) = orig.col(i);
                tmp(orig.n_rows) = nk;
                output.insert_cols(0,tmp);
            }
        }
    }
    else
    {
        for (Int nk=0; nk<=L; nk++)
        {
            IntVec tmp(1);
            tmp(0) = nk;
            output.insert_cols(0,tmp);
        }
    }
    return output;
}

IntMat get_indices(const Int & K, const Int & L)
{
    IntMat output(0,0,a::fill::zeros);
    for (int k=0; k<=K; k++)
    {
        output = return_larger_indices(output,L);
    }
    return output;
}

Int SysBase::z(const IntVec & intvec)
{
    int n{};
    int result{0};
    n = a::sum(intvec);
    //std::cout << " shape of flattening " << flattening.n_rows << " " << flattening.n_cols << "\n";
    for (int k=cK; 0<=k; k--)
    {
    //std::cout << " accessing " << n << " " << k << "\n";
        result += flattening(n,k);
        n -= intvec(k);
    }
    return result;
}


SysBase::SysBase(const Input &inp)
{
    // Getting MPI variables
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    bath_type = inp.bath_type;
    heom_switch = inp.heom_switch;
    t0 = inp.t0;
    t1 = inp.t1;
    dt = inp.dt;
    dt_sample = inp.dt_sample;
    sample_step = std::lrint(dt_sample/dt);
    if (sample_step == 0) sample_step=1;
    dt_sample = sample_step*dt;
    n_t = static_cast<Int>((t1-t0)/dt);

    tau = static_cast<Int>(dt_sample/dt);
    if (tau==0) tau = 1;

    ts = Vec(n_t+1,a::fill::zeros);
    if (t0<0 && t1>0)
    {
        Int n_neg{static_cast<Int>(-t0/dt)};
        Int cnt{0};
        for (Int i=-n_neg; i<=n_t-n_neg; i++)
        {
            ts(cnt) = dt*i;
            cnt++;
        }
        t_offset = n_neg%tau;
    }
    else
    {
        for (Int i=0; i<=n_t; i++)
        {
            ts(i) = dt*i + t0;
        }
        t_offset = 0;
    }

    beta = inp.beta;
    mass = inp.mass;
    eta = inp.eta;
    omega_c = inp.omega_c;
    cK = inp.cK;
    cL = inp.cL;
    n_ee = inp.n_ee;
    n_ados = binom(cL+cK+1,cK+1);
    indices = get_indices(cK,cL);
    
    mpi_split = IntMat(2,n_procs,a::fill::zeros);
    // This gives the first and the last index for each thread
    int per_proc{n_ados/n_procs};
    int remainder{n_ados%n_procs};
    mpi_split(0,0) = 0;
    if (remainder>0)
    {
        mpi_split(1,0) = per_proc;
    }
    else
    {
        mpi_split(1,0) = per_proc-1;
    }
    for (int i=1; i<n_procs; i++)
    {
        mpi_split(0,i) = mpi_split(1,i-1)+1;
        if (i<remainder)
        {
            mpi_split(1,i) = mpi_split(0,i)+per_proc;
        }
        else
        {
            mpi_split(1,i) = mpi_split(0,i)+per_proc-1;
        }
    }
    //std::cout << "n_ados " << n_ados << "\n";
    //std::cout << "mpi_split\n";
    //std::cout << mpi_split << "\n";
    //std::exit(6);
    //
    //std::cout << "n_ados " << n_ados << "\n";
    //std::cout << "n_cols " << indices.n_cols << "\n";
    //std::cout << "indices\n";
    //for (int i=0; i<indices.n_cols; i++)
    //{
    //    std::cout << (indices.col(i)).t() << "\n";
    //}

    // Externally loaded variables
    cs.load("tmp_cs",a::raw_ascii);
    gammas.load("tmp_gammas",a::raw_ascii);
    rho0.load("tmp_rho0",a::raw_ascii);
    ham.load("tmp_ham",a::raw_ascii);
    white_ham.load("tmp_white_ham",a::raw_ascii);
    q_mat.load("tmp_q_mat",a::raw_ascii);
    qq_mat = q_mat*q_mat;
    dq_mat.load("tmp_dq_mat",a::raw_ascii);
    flattening.load("tmp_flattening",a::raw_ascii);
    CplxVec low_T_coefs;
    low_T_coefs.load("tmp_low_T_coefs",a::raw_ascii);
    low_T_coef_R = low_T_coefs(0);
    low_T_coef_I = low_T_coefs(1);

    iham = (-imu/hbar)*ham;
    iq_mat = (-imu/hbar)*q_mat;

    ados = CplxCube(n_ee,n_ee,n_ados,a::fill::zeros);
    k1 = CplxCube(n_ee,n_ee,n_ados,a::fill::zeros);
    k2 = CplxCube(n_ee,n_ee,n_ados,a::fill::zeros);
    k3 = CplxCube(n_ee,n_ee,n_ados,a::fill::zeros);
    k4 = CplxCube(n_ee,n_ee,n_ados,a::fill::zeros);
    ados.slice(0) = rho0;

    if (rank==0)
    {
        std::cout << "Contents of SysBase:\n";
        std::cout << "heom_switch         " << heom_switch << '\n';
        std::cout << "bath_type           " << bath_type << '\n';
        std::cout << "\n";
        std::cout << "Simulation parameters\n";
        std::cout << "t0                  " << t0 << '\n';
        std::cout << "t1                  " << t1 << '\n';
        std::cout << "dt                  " << dt << '\n';
        std::cout << "dt_sample           " << dt_sample << '\n';
        std::cout << "\n";
        std::cout << "System parameters\n";
        std::cout << "beta                " << beta << '\n';
        std::cout << "mass                " << mass << '\n';
        std::cout << "eta                 " << eta << '\n';
        std::cout << "omega_c             " << omega_c << '\n';
        std::cout << "cK                  " << cK << '\n';
        std::cout << "cL                  " << cL << '\n';
        std::cout << "\n";
        std::cout << "\n";
    }
}
//****************************************************************************

void SysBase::change_to_white_ham()
{
    tmp_ham = ham;
    ham = white_ham;
    iham = (-imu/hbar)*ham;
}
void SysBase::change_from_white_ham()
{
    ham = tmp_ham;
    iham = (-imu/hbar)*ham;
}
Real SysBase::get_magnitude()
{
    tmpcplx = a::trace(ados.slice(0)*q_mat*q_mat);
    return std::abs(tmpcplx);
}
void SysBase::write_tcf_qq(const String & filename, const Real & time)
{
    tmpcplx = a::trace(ados.slice(0)*q_mat);
    write_out(filename+"_real_au_dat", time, std::real(tmpcplx) );
    write_out(filename+"_imag_au_dat", time, std::imag(tmpcplx) );
}
void SysBase::write_tcf_q2q2(const String & filename, const Real & time)
{
    tmpcplx = a::trace(ados.slice(0)*q_mat*q_mat);
    write_out(filename+"_real_au_dat", time, std::real(tmpcplx) );
    write_out(filename+"_imag_au_dat", time, std::imag(tmpcplx) );
}
void SysBase::write_tcf_q4q4(const String & filename, const Real & time)
{
    tmpcplx = a::trace(ados.slice(0)*q_mat*q_mat*q_mat*q_mat);
    write_out(filename+"_real_au_dat", time, std::real(tmpcplx) );
    write_out(filename+"_imag_au_dat", time, std::imag(tmpcplx) );
}
bool SysBase::check_stability(const Real & maxval)
{
    stability = true;
    for (int i=0; i<n_ee; i++)
    {
        if (std::abs(ados.slice(0)(i,i))>maxval)
        {
            stability = false;
        }
    }
    return stability;
}
void SysBase::postmultiply_by_qs()
{
    CplxMat tmp;
    tmp.copy_size(q_mat);
    tmp.fill(0.0);
    tmp.set_real(q_mat);
    ados.each_slice() *= tmp;
}
void SysBase::propagate()
{
    //std::cout << "*************************************************************\n";
    //std::cout << "*************************************************************\n";
    derivative(ados, dt/2.0, k1);
    //std::cout << "k1.slice(0)\n" << k1.slice(0) << "\n";
    derivative(ados+k1, dt/2.0, k2);
    //std::cout << "k2.slice(0)\n" << k2.slice(0) << "\n";
    derivative(ados+k2, dt, k3);
    //std::cout << "k3.slice(0)\n" << k3.slice(0) << "\n";
    derivative(ados+k3, dt/6.0, k4);
    //std::cout << "k4.slice(0)\n" << k4.slice(0) << "\n";
    ados += k1/3.0 + k2*(2.0/3.0) + k3/3.0 + k4;
    //std::cout << "ados.slice(0)\n" << ados.slice(0) << "\n";
    //std::cout << "*************************************************************\n";
    //std::cout << "*************************************************************\n";
}

void SysBase::derivative(const CplxCube & ado, const Real & coef, CplxCube & drdt)
{
    drdt.fill(0.0);
    for (int i=mpi_split(0,rank); i<=mpi_split(1,rank); i++)
    {
        curvec = indices.col(i);
        curind = z(curvec);
        rho = ado.slice(curind)*coef;
        if (bath_type == NO_BATH)
        {
            drdt.slice(curind) += iham*rho - rho*iham;
            // NB: iham includes -i/hbar
        }
        else
        {
            //if elif bath_type in {"debye","debye_correction","debye_correction2","white"}:
            if (bath_type == DEBYE)
            {
                // Terms: Liouvillian + friction i.e. n = n
                drdt.slice(curind) += iham*rho - rho*iham
                                    - a::dot(curvec, gammas)*rho
                                        ;
            }
            else if (bath_type == DEBYE_CORRECTION)
            {
                drdt.slice(curind) += iham*rho - rho*iham
                                    - a::dot(curvec, gammas)*rho
                                    - low_T_coef_R*(
                                            qq_mat*rho
                                            - 2*q_mat*rho*q_mat
                                            + rho*qq_mat)
                                        ;
            }
            else if (bath_type == DEBYE_CORRECTION2)
            {
                drdt.slice(curind) += iham*rho - rho*iham
                                    - a::dot(curvec, gammas)*rho
                                    - (low_T_coef_R+low_T_coef_I)*qq_mat*rho
                                    - (low_T_coef_R-low_T_coef_I)*rho*qq_mat
                                    + 2.0*low_T_coef_R*q_mat*rho*q_mat
                                        ;
            }
            else if (bath_type == WHITE)
            {
                drdt.slice(curind) += iham*rho - rho*iham
                                    - a::dot(curvec, gammas)*rho
                                    - low_T_coef_R*(
                                            qq_mat*rho
                                            - 2*q_mat*rho*q_mat
                                            + rho*qq_mat
                                            )
                                    - low_T_coef_I*(
                                            q_mat*dq_mat*rho
                                            +q_mat*rho*dq_mat
                                            -dq_mat*rho*q_mat
                                            -rho*dq_mat*q_mat
                                            //q_mat*dq_mat*rho
                                            //-q_mat*rho*dq_mat
                                            //-dq_mat*rho*q_mat
                                            //+rho*dq_mat*q_mat
                                            )
                                        ;
            }
            // Calculating the "Commutator and anticommutator"
            // NB: iqs = -i/hbar * qs
            q_rho = iq_mat*rho;
            rho_q = rho*iq_mat;
            comm = q_rho-rho_q;
            // rho_+k
            if (a::sum(curvec)>0)
            {
                for (int k=0; k<=cK; k++)
                {
                    if (curvec(k)!=0)
                    {
                        nextvec = curvec;
                        nextvec(k) -= 1;
                        drdt.slice(z(nextvec)) += comm;
                    }
                }
            }
            // rho_-k
            if (a::sum(curvec)<cL)
            {
                // k > 0
                for (int k=1; k<=cK; k++)
                {
                    nextvec = curvec;
                    nextvec(k) += 1;
                    drdt.slice(z(nextvec)) += static_cast<double>(nextvec(k))*cs(k)*comm;
                }
                // k = 0
                nextvec = curvec;
                nextvec(0) += 1;
                drdt.slice(z(nextvec)) += nextvec(0)*(
                        cs(0)*q_rho-std::conj(cs(0))*rho_q);
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, drdt.memptr(), drdt.n_elem, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
}
