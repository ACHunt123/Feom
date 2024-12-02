// File: input.cpp
#include "general.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include "input.h"
#include "utils.h"

template <typename T>
using maptype = const std::map<String,T>;
template <typename T>
T value_from_map(String &key, maptype<T> &mymap)
{
    if (mymap.count(key)>0) return mymap.at(key);
    else
    {
        std::cout << "Unable to find keyword: " << key << "\n";
        std::exit(1);
    } 
}
template <typename T>
String key_from_map(T &inpval, maptype<T> &mymap)
{
    for (const auto& [key, value] : mymap)
    {
        if (value == inpval) return key;
    }
    std::cout << "Unable to find value: " << inpval << "\n";
    std::exit(1);
    return "";
}

Input::Input(String filename)
{
    std::ifstream input_file (filename);
    if (!input_file.is_open())
    {
        std::cout << "Unable to open the input file";
        std::exit(1);
    }
    else
    {
        String current_line{};
        String type_temp;
        // Read the file
        //
        //
        //
        //
        getline(input_file, current_line);//  heompy input file
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  ***Checkpoints***
        getline(input_file, current_line);//  load_checkpoint               
        getline(input_file, current_line);//  save_checkpoint               
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  ***TO be calculated***
        getline(input_file, current_line);//  heom_switch                   
        extract_line(current_line, type_temp);
        heom_switch = value_from_map(type_temp,heom_switch_map);
        getline(input_file, current_line);//  traj_switch                   
        getline(input_file, current_line);//  ado_save                      
        getline(input_file, current_line);//  ado_trend_save                
        getline(input_file, current_line);//  write_density                 
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  ***Simlation switches***
        getline(input_file, current_line);//  sim_type                      
        getline(input_file, current_line);//  bath_type                     
        extract_line(current_line, type_temp);
        bath_type = value_from_map(type_temp,bath_type_map);
        getline(input_file, current_line);//  sop_decomposition_type        
        getline(input_file, current_line);//  init_state                    
        getline(input_file, current_line);//  pot_type                      
        getline(input_file, current_line);//  a_ren_type                    
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  **Numeric parameters**
        getline(input_file, current_line);//  stability_threshold           
        extract_line(current_line, stability_threshold);
        getline(input_file, current_line);//  truncation_type               
        getline(input_file, current_line);//  scaling_type                  
        getline(input_file, current_line);//  pruning_type                  
        getline(input_file, current_line);//  cL                            
        extract_line(current_line, cL);
        getline(input_file, current_line);//  cK                            
        extract_line(current_line, cK);
        getline(input_file, current_line);//  cutoff_coef                   
        getline(input_file, current_line);//  n_EE                          
        extract_line(current_line, n_ee);
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  **System parameters**
        getline(input_file, current_line);//  mass                          
        extract_line(current_line, mass);
        getline(input_file, current_line);//  beta                          
        extract_line(current_line, beta);
        getline(input_file, current_line);//  eta                           
        getline(input_file, current_line);//  gamma                         
        extract_line(current_line, omega_c);
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  **Simulation parameters**
        getline(input_file, current_line);//  q0                            
        getline(input_file, current_line);//  q1                            
        getline(input_file, current_line);//  dq                            
        getline(input_file, current_line);//  t0                            
        extract_line(current_line, t0);
        getline(input_file, current_line);//  t1                            
        extract_line(current_line, t1);
        getline(input_file, current_line);//  dt                            
        extract_line(current_line, dt);
        getline(input_file, current_line);//  t_sample                      
        extract_line(current_line, dt_sample);
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  ***Wavepacket parameters***
        getline(input_file, current_line);//  wp_q0                         
        getline(input_file, current_line);//  wp_sigma                      
        getline(input_file, current_line);//  wp_p0                         
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  ***Trajectory parameters***
        getline(input_file, current_line);//  rand_num_seed                 
        getline(input_file, current_line);//  n_trajs                       
        getline(input_file, current_line);//  tr_dt                         
        getline(input_file, current_line);//  n_b                           
        getline(input_file, current_line);//  n_therm                       
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  ***Parameters of potentials***
        getline(input_file, current_line);//  ho_omega                      
        getline(input_file, current_line);//  ho_x0                         
        getline(input_file, current_line);//  ho_x0_switch                  
        getline(input_file, current_line);//  ho_const                      
        getline(input_file, current_line);//  ho_const_switch               
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  tani_omega_a                  
        getline(input_file, current_line);//  tani_f_coef                   
        getline(input_file, current_line);//  tani_f_coef_switch            
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  morse_k                       
        getline(input_file, current_line);//  morse_de                      
        getline(input_file, current_line);//  morse_x0                      
        getline(input_file, current_line);//  morse_x0_switch               
        getline(input_file, current_line);//  morse_const                   
        getline(input_file, current_line);//  morse_const_switch            
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  poly_coefs                    
        getline(input_file, current_line);//  poly_coefs_switch             
        getline(input_file, current_line);//  poly_x0                       
        getline(input_file, current_line);//  poly_x0_switch                
        getline(input_file, current_line);//  
        getline(input_file, current_line);//  anh_c                         
        //extract_line(current_line, random_seed);

        //getline(input_file, current_line);
        //extract_line(current_line, n_trajs);

        //getline(input_file, current_line);
        //extract_line(current_line, n_trajs_iter);

        //getline(input_file, current_line);
        //extract_line(current_line, max_n_iter);

        //getline(input_file, current_line);
        //extract_line(current_line, error_threshold);

        //getline(input_file, current_line);
        //extract_line(current_line, therm_time);

        //getline(input_file, current_line);
        //extract_line(current_line, n_therm_cycles);

        //getline(input_file, current_line);
        //extract_line(current_line, n_relax_cycles);

        //getline(input_file, current_line);
        //getline(input_file, current_line);
        //// Simulation parameters
        //getline(input_file, current_line);
        //extract_line(current_line, t0);

        //getline(input_file, current_line);
        //extract_line(current_line, t1);

        //getline(input_file, current_line);
        //extract_line(current_line, dt);

        //getline(input_file, current_line);
        //extract_line(current_line, dt_sample);

        //getline(input_file, current_line);
        //getline(input_file, current_line);
        //// System parameters
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //sim_type = value_from_map(type_temp,sim_type_map);

        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //system_type = value_from_map(type_temp,system_type_map);

        //getline(input_file, current_line);
        //extract_line(current_line, beta);

        //getline(input_file, current_line);

        //getline(input_file, current_line);
        //getline(input_file, current_line);
        //// Bath parameters
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //bath_type = value_from_map(type_temp,bath_type_map);

        //getline(input_file, current_line);
        //extract_line(current_line, eta);
        //// For all current definitions of the baths zeta0=eta
        //zeta0 = eta;

        //getline(input_file, current_line);
        //extract_line(current_line, omega_c);

        //getline(input_file, current_line);
        //getline(input_file, current_line);
        //// Paremeters of SysTrajs/SysBerk
        //getline(input_file, current_line);
        //extract_line(current_line, n_bath_modes);
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //noise_type = value_from_map(type_temp,noise_type_map);
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //distr_type = value_from_map(type_temp,distr_type_map);

        //getline(input_file, current_line);
        //getline(input_file, current_line);
        //// Paremeters of SysMats
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //matspot_type = value_from_map(type_temp,matspot_type_map);
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //mf_type = value_from_map(type_temp,mf_type_map);
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //ac_type = value_from_map(type_temp,ac_type_map);
        //getline(input_file, current_line);
        //extract_line(current_line, n_beads);
        //getline(input_file, current_line);
        //extract_line(current_line, n_mats_modes);
        //getline(input_file, current_line);
        //extract_line(current_line, adiabatic_separation);
        //getline(input_file, current_line);
        //extract_line(current_line, sdt);

        //getline(input_file, current_line);
        //getline(input_file, current_line);
        //// Parameters of potentials
        //getline(input_file, current_line);
        //extract_line(current_line, type_temp);
        //pot_type = value_from_map(type_temp,pot_type_map);

        //getline(input_file, current_line);
        //extract_line(current_line, ho_omega);

        //getline(input_file, current_line);
        //extract_line(current_line, ho_q0);

        //getline(input_file, current_line);
        //extract_line(current_line, ho_const);

        //getline(input_file, current_line);
        //extract_line(current_line, hydoh_c);
    }
}

void Input::write_log(String logfile)
{
    std::ofstream log(logfile,std::ios::app);
    if (!log.is_open())
    {
        // Print an error and exit
        std::cerr << "Uh oh, "<< logfile <<" could not be opened for writing!"
                                                                << std::endl;
        std::exit(1);
    }
    output_data(log);
    log.close();
}

void Input::output_data(std::ostream & log)
{
    log << std::setprecision(17);
    log << "Logfile of a heomc calculation\n";
    log << "\n";
    //log << "***Numeric parameters***\n";
    //log << "random_seed                      " << random_seed << '\n';
    //log << "n_trajs                          " << n_trajs << '\n';
    //log << "n_trajs_iter                     " << n_trajs_iter << '\n';
    //log << "max_n_iter                       " << max_n_iter << '\n';
    //log << "error_threshold                  " << error_threshold << '\n';
    //log << "therm_time/(beta*hbar)           " << therm_time << '\n';
    //log << "n_therm_cycles                   " << n_therm_cycles << '\n';
    //log << "n_relax_cycles                   " << n_relax_cycles << '\n';
    //log << "\n";
    //log << "***Simulation parameters***\n";
    log << "t0                               " << t0 << '\n';
    log << "t1                               " << t1 << '\n';
    log << "dt                               " << dt << '\n';
    log << "dt_sample                        " << dt_sample << '\n';
    log << "\n";
    //log << "***System parameters***\n";
    //log << "sim_type                         "
    //    << key_from_map(sim_type,sim_type_map) << '\n';
    //log << "system_type                      "
    //    << key_from_map(system_type,system_type_map) << '\n';
    log << "beta                             " << beta << '\n';
    log << "mass                             " << mass << '\n';
    log << "\n";
    //log << "***Bath parameters***\n";
    log << "heom_switch                      "
        << key_from_map(heom_switch,heom_switch_map) << '\n';
    log << "bath_type                        "
        << key_from_map(bath_type,bath_type_map) << '\n';
    log << "eta                              " << eta << '\n';
    log << "omega_c                          " << omega_c << '\n';
    log << "cK                               " << cK << '\n';
    log << "cL                               " << cL << '\n';
    log << "n_ee                             " << n_ee << '\n';
    //log << "\n";
    //log << "***SysTrajs/SysBerk parameters***\n";
    //log << "n_bath_modes                     " << n_bath_modes << '\n';
    //log << "noise_type                       "
    //    << key_from_map(noise_type,noise_type_map) << '\n';
    //log << "distr_type                       "
    //    << key_from_map(distr_type,distr_type_map) << '\n';
    //log << "\n";
    //log << "***SysRing parameters***\n";
    //log << "matspot_type                        "
    //    << key_from_map(matspot_type,matspot_type_map) << '\n';
    //log << "n_beads                          " << n_beads << '\n';
    //log << "n_mats_modes                     " << n_mats_modes << '\n';
    //log << "adiabatic_separation             " << adiabatic_separation << '\n';
    //log << "sdt                             " << sdt << '\n';
    //log << "\n";
    //log << "***Parameters of potentials***\n";
    //log << "pot(ential)_type                 "
    //    << key_from_map(pot_type,pot_type_map) << '\n';
    //log << "ho_omega                         " << ho_omega << '\n';
    //log << "ho_q0                            " << ho_q0 << '\n';
    //log << "ho_const                         " << ho_const << '\n';
    //log << "hydoh_c                          " << hydoh_c << '\n';
    //log << "\n";
    //log << "\n";
    //log << "Created in Input object\n";
    //log << "zeta0                            " << zeta0 << '\n';
    log << "\n";
}
