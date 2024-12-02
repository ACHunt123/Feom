#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
import pickle

import heom.dvr as dvr
from heom.input import InputObj
from heom.type_dicts import potentials, states
from heom.general import hbar, pi, formatflt, write_complex_pair


def main():
    file_input = sys.argv[1]
    if len(sys.argv)>2:
        name = "_" + sys.argv[2]
    else:
        name = ""
    file_tcf_qq_real = "tcf_qq"+name+"_stnd_real_au_dat"
    file_tcf_qq_imag = "tcf_qq"+name+"_stnd_imag_au_dat"
    file_kubo_tcf_qq_real = "tcf_qq"+name+"_kubo_real_au_dat"
    file_kubo_tcf_qq_imag = "tcf_qq"+name+"_kubo_imag_au_dat"
    file_tcf_q2q2_real = "tcf_q2q2"+name+"_stnd_real_au_dat"
    file_tcf_q2q2_imag = "tcf_q2q2"+name+"_stnd_imag_au_dat"
    file_kubo_tcf_q2q2_real = "tcf_q2q2"+name+"_kubo_real_au_dat"
    file_kubo_tcf_q2q2_imag = "tcf_q2q2"+name+"_kubo_imag_au_dat"
    inp = InputObj(file_input)
    pot = potentials[inp.pot_type](inp)

    #****************************************************************************
    # If renormalisation potential is to be included
    #pot = potentials["PotEff"](pot,a_ren=inp.eta*inp.gamma)
    #****************************************************************************

    state = states[inp.init_state](inp,pot)

    q1 = inp.q1
    q0 = inp.q0
    dq = inp.dq
    mass = inp.mass
    beta = inp.beta
    t0 = inp.t0
    t1 = inp.t1
    t_sample = inp.t_sample
    n_EE = inp.n_EE
    print("Using {} lowest energy wave functions".format(n_EE))

    n_q = int((q1-q0)/dq)
    qs = np.empty([n_q])
    for i in range(n_q):
        # + 1 is there just for consistency with Fortran
        qs[i] = dq*(i + 1 + int(q0/dq))

    #kinE = -(hbar**2/(2*mass))*dvr.derivative2(qs,dq)
    #potE = np.identity(n_q)
    #for i in range(n_q):
    #    potE[i,i] = pot.calc(qs[i])
    #ham = kinE + potE
    ham = dvr.hamiltonian(qs, pot, mass)

    evals, evecs = np.linalg.eigh(ham)
    # Normalising Wavefunctions
    for n in range(len(evals)):
        evecs[:,n] /= np.sqrt(dq*np.vdot(evecs[:,n],evecs[:,n])) 

    if len(evals)<n_EE:
        n_EE = len(evals)

    print("Number of used functions", n_EE)
    ## Plotting wavefunctions
    #fig = plt.figure(figsize=(8, 6), dpi=160)
    #for n in range(n_EE):
    #    plt.plot(qs,np.real(evecs[:,n]), label=str(n)+" real")
    #    #plt.plot(qs,np.imag(evecs[:,n]), label=str(n)+" imaginary", ls="--")
    #plt.plot(qs, pot.calc(qs), label="Potential", c='k')
    #plt.legend()
    #plt.show()
    #plt.clf()
    #plt.close()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('position')
    ax1.set_ylabel('wave function', color='C0')
    for n in range(n_EE):
        ax1.plot(qs,np.real(evecs[:,n]), label=str(n)+" real")
        #plt.plot(qs,np.imag(evecs[:,n]), label=str(n)+" imaginary", ls="--")
    ax1.tick_params(axis='y', labelcolor='C0')
    plt.legend()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'k'
    ax2.set_ylabel('potential', color=color)  # we already handled the x-label with ax1
    ax2.plot(qs, pot.calc(qs), label="Potential", c=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

    # Calculating the Q_mn matrix elements
    q_mat = np.zeros([n_EE,n_EE])
    for n in range(n_EE):
        for m in range(n_EE):
            q_mat[n,m] = dq*np.vdot(evecs[:,n],qs*evecs[:,m])
    # Calculating the Q2_mn matrix elements
    q2_mat = np.zeros([n_EE,n_EE])
    for n in range(n_EE):
        for m in range(n_EE):
            q2_mat[n,m] = dq*np.vdot(evecs[:,n],np.square(qs)*evecs[:,m])
    # Calculating the partition function
    partfun = 0
    for n in range(n_EE):
        partfun += np.exp(-beta*evals[n])
    # Calculating TCF
    def make_tcf(t,mat):
        cf = 0
        for n in range(n_EE):
            for m in range(n_EE):
                cf += mat[n,m]*mat[m,n]*np.exp(
                        -beta*evals[n]-(1j*t/hbar)*(evals[n]-evals[m]))
        return cf/partfun

    def make_kubo_tcf(t,mat):
        cf = 0
        for n in range(n_EE):
            for m in range(n_EE):
                tmp = (mat[n,m]*mat[m,n]
                    *np.exp(-(1j*t/hbar)*(evals[n]-evals[m])))
                if m != n:
                    tmp *= ((np.exp(-beta*evals[m])-np.exp(-beta*(evals[n])))
                            / (evals[n]-evals[m]))
                elif m == n:
                    tmp *= beta*np.exp(-beta*evals[m])
                cf += tmp

        return cf/(partfun*beta)

    ts = np.linspace(0,t1, int(t1/t_sample))
    tcf = []
    for i in range(len(ts)):
        message = ("Making TCF 〈 qq(t) 〉 step {:8d} out of {:8d} steps"
                    ).format(i,len(ts))
        sys.stdout.write("\r"+message)
        tmp = make_tcf(ts[i],q_mat)
        tcf.append(tmp)
        write_complex_pair(file_tcf_qq_real,file_tcf_qq_imag,ts[i],tmp)
    print("")

    kubo_tcf = []
    for i in range(len(ts)):
        message = ("Making Kubo TCF 〈 qq(t) 〉 step {:8d} out of {:8d} steps"
                    ).format(i,len(ts))
        sys.stdout.write("\r"+message)
        tmp = make_kubo_tcf(ts[i],q_mat)
        kubo_tcf.append(tmp)
        write_complex_pair(file_kubo_tcf_qq_real,file_kubo_tcf_qq_imag,ts[i],tmp)
    print("")
    for i in range(len(ts)):
        message = ("Making TCF 〈 q²q²(t) 〉 step {:8d} out of {:8d} steps"
                    ).format(i,len(ts))
        sys.stdout.write("\r"+message)
        tmp = make_tcf(ts[i],q2_mat)
        tcf.append(tmp)
        write_complex_pair(file_tcf_q2q2_real,file_tcf_q2q2_imag,ts[i],tmp)
    print("")

    kubo_tcf = []
    for i in range(len(ts)):
        message = ("Making Kubo TCF 〈 q²q²(t) 〉 step {:8d} out of {:8d} steps"
                    ).format(i,len(ts))
        sys.stdout.write("\r"+message)
        tmp = make_kubo_tcf(ts[i],q2_mat)
        kubo_tcf.append(tmp)
        write_complex_pair(file_kubo_tcf_q2q2_real,file_kubo_tcf_q2q2_imag,ts[i],tmp)
    print("")
    quit()

    # Plotting TCF
    #fig = plt.figure(figsize=(8, 6), dpi=160)
    #plt.plot(ts,tcf)
    #plt.title("QQ(t) TCF")
    #plt.show()
    #plt.clf()
    #plt.close()

    def create_rho0():
        rho = np.zeros([n_q,n_q],dtype=complex)
        for n in range(n_EE):
            rho += (np.exp(-beta*evals[n])/partfun)*np.outer(evecs[:,n],evecs[:,n])
        return rho

    # Exporting equilibrium density matrix
    rho0 = create_rho0()
    with open("rho0", 'wb') as f:
        pickle.dump((qs,rho0),f)
    


    def check_init_state(rho):
        print("Checking if the initial state is the equilibrium state for the potential")
        #rho = np.real(rho)
        #rho /= np.trace(rho)
        rho0 = state.rho(qs)
        #rho0 = np.real(rho0)
        #rho0 /= np.trace(rho0)

        good = 0
        bad = 0
        tolerance = 1e-12
        for i in range(n_EE):
            for j in range(n_EE):
                if abs(rho[i,j]-rho0[i,j]) < tolerance*abs(rho0[i,j]):
                    good += 1
                else:
                    bad += 1
        print("Out of {} points, {} %, within {} % of each other.".format(good+bad, 100*good/(good+bad), 100*tolerance))
        fig = plt.figure(figsize=(6, 6), dpi=160)
        plt.contour(qs, qs, rho, cmap="Reds", linewidths=3)
        plt.contour(qs, qs, rho0, cmap="Greens", linestyles='dashed', linewidths=1)
        plt.show()

    check_init_state(rho0)

if __name__=="__main__":
    main()
    print("Program finished succesfully")
