// File: input.h
#ifndef INPUT_H
#define INPUT_H
#include "general.h"

class Input
{
public:
    Real t0{};
    Real t1{};
    Real dt{};
    Real dt_sample{};
    Real stability_threshold{};
    Bath_type bath_type;
    Heom_switch heom_switch;
    Real eta{};
    Real omega_c{};
    Real beta{};
    Real mass{};
    Int cL{};
    Int cK{};
    Int n_ee{};

//    // Numeric parameters
//    Int random_seed{};
//    Int n_trajs{};
//    Int n_trajs_iter{};
//    Int max_n_iter{};
//    Real error_threshold{};
//    Real therm_time{};
//    Real n_therm_cycles{};
//    Real n_relax_cycles{};
//    // Simulation parameters
//    // System parameters
//    Sim_type sim_type;
//    System_type system_type;
//    Noise_type noise_type;
//    Distr_type distr_type;
//    // Bath parameters
//    Real zeta0{};
//    // SysTrajs/SysBerk parameters
//    Int n_bath_modes{};
//    // SysMats
//    Matspot_type matspot_type;
//    MF_type mf_type;
//    AC_type ac_type;
//    Int n_beads{};
//    Int n_mats_modes{};
//    Real adiabatic_separation{};
//    Real sdt{};
//    // Parameters of the potentials
//    Pot_type pot_type;
//    Real ho_omega{};
//    Real ho_q0{};
//    Real ho_const{};
//    Real hydoh_c{};
    // Member functions
    Input(String filename);
    void write_log(String logfile);
    void output_data(std::ostream & log);
};
#endif
