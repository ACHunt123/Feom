// File: heomc.cpp
#include "general.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
// For signal handling:
#include <csignal>
// My files
#include "general.h"
#include "systems.h"
#include "input.h"
#include "timer.h"
#include "utils.h"

String file_tcf_qq{"tcf_qq_"};
String file_tcf_q2q2{"tcf_q2q2_"};
String file_equilibration{"eq_"};
String file_q2{"q2_"};
String file_q4{"q4_"};

void write_out(SysBase & sys, const Int & t, const String & name)
{
    if (sys.ts(t)<0)
    {
        sys.write_tcf_q2q2(file_q2+name, sys.ts(t));
        sys.write_tcf_q4q4(file_q4+name, sys.ts(t));
    }
    if (sys.heom_switch == HEOM_TCF_QQ or sys.heom_switch == HEOM_KUBO_QQ)
    {
        if (sys.ts(t)<=0)
        {
            sys.write_tcf_q2q2(file_equilibration+name, sys.ts(t));
        }
        else
        {
            sys.write_tcf_qq(file_equilibration+name, sys.ts(t));
        }
    }
    else //i.e. calculating quadratic TCF
    {
        if (sys.ts(t)<=0)
        {
            sys.write_tcf_q4q4(file_equilibration+name, sys.ts(t));
        }
        else
        {
            sys.write_tcf_q2q2(file_equilibration+name, sys.ts(t));
        }
    }
    //if inp.write_density!=0:
    //    sys.write_density(file_density, sys.ts(t))
    //    if type(sys) is sysulations["lowT_qm_ee"]:
    //        sys.write_density_EE(file_density_EE, sys.ts(t))
    //if inp.ado_trend_save == 1:
    //    sys.write_ado_trend(file_ado_trend, sys.ts(t))
}

int main (int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank{0};
    int n_procs{0};
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    signal(SIGINT, signal2exitfile);  
    signal(SIGTERM, signal2exitfile);  
    Timer timer;
    timer.start();
    if (rank==0) std::cout << "Running program heomc by Adam Prada.\n";
    String file_input{};
    String name{"heomc"};
    String logfile{"log"};
    if (argc < 2)
    {
        std::cerr << "Please give the input file name as an argument.\n";
        std::exit(1);
    }
    else file_input = argv[1];
    if (argc > 2)
    {
        name = argv[2];
        if (rank==0) std::cout << "The detected calculation name is: " << name << "\n";
    }
    // Initialising objects
    Input inp{file_input};
    if (rank==0) inp.write_log(logfile);
    SysBase sys(inp);

    bool delta_kick{true};

    bool dirprod{true};
    if (sys.ts(0)<0) dirprod = false;
    if (rank==0) std::cout << "dirprod: " << dirprod << "\n";
    bool stability{true}, exitbool{false};
    Int t{0},t_zeropoint{0};
    //***************************************************************************
    // Propagating for t<0
    //***************************************************************************
    
    if (delta_kick and not dirprod and inp.bath_type==WHITE)
    {
        sys.change_to_white_ham();
    }
    while (sys.ts(t)<0 and t<sys.n_t && stability and not exitbool)
    {
        sys.propagate();
        t++;
        if (delta_kick and not dirprod and inp.bath_type==WHITE and t == t_zeropoint+1)
        {
            sys.change_from_white_ham();
        }
        stability = check_stability(inp.stability_threshold,sys.get_magnitude());
        exitbool = check_exit();
        if (rank==0 and t%sys.tau==sys.t_offset) write_out(sys,t,name);
    }
    //***************************************************************************
    // Making change at t=0
    //***************************************************************************
    if (sys.heom_switch==HEOM_TCF_QQ)
    {
        // <qq(t)>
        sys.postmultiply_by_qs();
        // Writing out t=0 tcf
        if (rank==0 and std::abs(sys.ts(t))<1e-7)
        {
            sys.write_tcf_qq(file_tcf_qq+name+"_stnd", sys.ts(t));
        }
    }
    else if (sys.heom_switch==HEOM_KUBO_QQ)
    {
        // Kubo <qq(t)>
        // Writing out t=0 tcf
        if (rank==0 and std::abs(sys.ts(t))<1e-7)
        {
            sys.write_tcf_qq(file_tcf_qq+name+"_kubo", sys.ts(t));
        }
    }
    else if (sys.heom_switch==HEOM_TCF_Q2Q2)
    {
        // <q^2q^2(t)>
        sys.postmultiply_by_qs();
        sys.postmultiply_by_qs();
        // Writing out t=0 tcf
        if (rank==0 and std::abs(sys.ts(t))<1e-7)
        {
            sys.write_tcf_q2q2(file_tcf_q2q2+name+"_stnd", sys.ts(t));
        }
    }
    else if (sys.heom_switch==HEOM_KUBO_Q2Q2)
    {
        // Kubo <q^2q^2(t)>
        // Writing out t=0 tcf
        if (rank==0 and std::abs(sys.ts(t))<1e-7)
        {
            sys.write_tcf_q2q2(file_tcf_q2q2+name+"_kubo", sys.ts(t));
        }
    }


    if (delta_kick and dirprod and inp.bath_type==WHITE)
    {
        sys.change_to_white_ham();
    }
    t_zeropoint = t;
    while (t<sys.n_t and stability and not exitbool)
    {
        sys.propagate();
        t++;
        if (delta_kick and dirprod and inp.bath_type==WHITE and t == t_zeropoint+1)
        {
            sys.change_from_white_ham();
        }
        stability = check_stability(inp.stability_threshold,sys.get_magnitude());
        exitbool = check_exit();
        if (rank==0 and t%sys.tau==sys.t_offset)
        {
            write_out(sys,t,name);
            if (sys.heom_switch==HEOM_TCF_QQ)
            {
                sys.write_tcf_qq(file_tcf_qq+name+"_stnd", sys.ts(t));
            }
            else if (sys.heom_switch==HEOM_KUBO_QQ)
            {
                sys.write_tcf_qq(file_tcf_qq+name+"_kubo", sys.ts(t));
            }
            else if (sys.heom_switch==HEOM_TCF_Q2Q2)
            {
                sys.write_tcf_q2q2(file_tcf_q2q2+name+"_stnd", sys.ts(t));
            }
            else if (sys.heom_switch==HEOM_KUBO_Q2Q2)
            {
                sys.write_tcf_q2q2(file_tcf_q2q2+name+"_kubo", sys.ts(t));
            }

        }
    }

    timer.end();

    if (rank==0)
    {
        std::cout << "\nTotal running time: " << timer.get_string() << "\n";
        write_out_default(logfile, String{"Total running time / s:"});
        write_out_default(logfile, timer.get_seconds());
        write_out_default("time_"+std::to_string(static_cast<Int>(timer.get_seconds())), timer.get_seconds());
        write_out_default(logfile, String{"\n"});
        write_out_default(logfile, String{"Total running time / human readable:"});
        write_out_default(logfile, timer.get_string());
    }
    MPI_Finalize();
}
