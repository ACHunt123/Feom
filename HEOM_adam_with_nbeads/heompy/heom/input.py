#!/usr/bin/env python
# File: input.py
"""Input object."""
import heom.type_dicts as tps

defaults = {
        "load_checkpoint":     0,
        "save_checkpoint":     1,
        "heom_switch":         "tcf_qq",
        "traj_switch":         0,
        "ado_save":            0,
        "ado_trend_save":      1,
        "write_density":       0,
        "sim_type":            "lowT_qm_ee",
        "bath_type":           "debye_correction",
        "sop_decomposition_type":"matsubara",
        "init_state":          "thermal",
        "pot_type":            "champagne",
        "a_ren_type":          "gamma*eta",
        "stability_threshold": 1000,
        "truncation_type":     "zero",
        "scaling_type":        "shi",
        "pruning_type":        "pruning",
        "cL":                  10,
        "cK":                  0,
        "cutoff_coef":         1.0e-5,
        "n_EE":                10,
        "mass":                1741.1, #1822.89
        "beta":                1052.5837645244185,
        "eta":                 18.2289,
        "gamma":               0.027337966,
        "q0":                  -0.7,
        "q1":                  0.8 ,
        "dq":                  0.05,
        "t0":                  0.0,
        "t1":                  862.8,
        "dt":                  1.0,
        "t_sample":            4.0,
        "wp_q0":               0.000,
        "wp_sigma":            0.700,
        "wp_p0":               0.000,
        "rand_num_seed":       1,
        "n_trajs":             5000,
        "tr_dt":               1.0e-3,
        "n_b":                 50,
        "n_therm":             20,
        "ho_omega":            0.013668983,
        "ho_x0":               0.000,
        "ho_x0_switch":        0.000,
        "ho_const":            0.000,
        "ho_const_switch":     0.000,
        "tani_omega_a":        2.000,
        "tani_f_coef":         0.000,
        "tani_f_coef_switch":  6.000,
        "morse_alpha":         4.000,
        "morse_de":            100.0,
        "morse_x0":            0.000,
        "morse_x0_switch":     1.500,
        "morse_const":         0.000,
        "morse_const_switch":  -4.50,
        "poly_coefs":          [0.000,0.000,2.2500,0.000],
        "poly_coefs_switch":   [0.000,-6.00,2.2500,0.000],
        "poly_x0":             0.000,
        "poly_x0_switch":      1.500,
        "anh_c":               2.000e-1,
}

defaults_rpsimc = {
                "random_seed": 485,
                    "n_trajs": 1000,
               "n_trajs_iter": 1000,
                         "dt": 0.1,
                        "sdt": 0.02,
                   "pot_type": "harmonic",
                 "max_n_iter": 1000,
            "error_threshold": 0.003,
     "therm_time/(beta*hbar)": 1.0,
             "n_therm_cycles": 20,
             "n_relax_cycles": 3,
     "pile_centroid_friction": 2*0.01703043563586307,
                         "t0": 0.0,
                         "t1": 1000.0,
                  "dt_sample": 1.0,
                   "sim_type": "frontier_modes",
                "system_type": "matsberk",
                       "beta": 2105.167529048837,
                       "mass": 1741.1,
                  "bath_type": "debye",
                        "eta": 59.30338297120238,
                    "omega_c": 2*0.01703043563586307,
               "n_bath_modes": 1000,
                 "noise_type": "real",
                 "distr_type": "dirprod",
               "matspot_type": "numeric",
                    "mf_type": "mf_on",
                    "ac_type": "ac_on",
                    "n_beads": 256,
               "n_mats_modes": 85,
       "adiabatic_separation": 16,
                   "ho_omega": 0.01703043563586307,
                      "ho_q0": 0.000,
                   "ho_const": 0.000,
                    "hydoh_c": 1.000,
}

class InputObj(object):
    def __init__(self, filename, checks=True):
        with open(filename, 'rt') as f:
            for li in f:
                line = li.split()
                if len(li)>1:
                    key = line[0]
                    val = line[1:]
                    if (key in defaults):
                        if type(defaults[key])==list:
                            setattr(self,key,[float(z) for z in val])
                        else:
                            setattr(self,key,type(defaults[key])(val[0]))
                    elif (key in defaults_rpsimc):
                        if type(defaults_rpsimc[key])==list:
                            setattr(self,key,[float(z) for z in val])
                        else:
                            setattr(self,key,type(defaults_rpsimc[key])(val[0]))
                    elif key[:2]=="**" or key=="rpsimc input file":
                        pass
                    else:
                        print("Ignoring key: ",key)
        # Discrepancy between rpsimc and heompy in name --- either one should be valid
        if hasattr(self,"dt_sample") and not hasattr(self,"t_sample"):
            setattr(self,"t_sample",self.dt_sample)
        for key in defaults:
            if not hasattr(self,key):
                setattr(self,key,defaults[key])

        # Running checks on input data (should implement)
        if checks:
            assert self.heom_switch in tps.heom_switches
            assert self.sim_type in tps.simulations
            assert self.bath_type in tps.bath_types
            assert self.init_state in tps.states
            assert self.pot_type in tps.potentials
            assert self.a_ren_type in tps.a_ren_types
            assert self.truncation_type in tps.truncation_types
            assert self.scaling_type in tps.scaling_types
            assert self.pruning_type in tps.pruning_types
            assert self.sop_decomposition_type in tps.sop_decomposition_types

if __name__=="__main__":
    import sys
    file_input = sys.argv[1]
    inp = InputObj(file_input)
    print("Object dictionary")
    print(inp.__dict__)
    print("Saved values")
    for key,val in inp.__dict__.items():
        print(key,val)
