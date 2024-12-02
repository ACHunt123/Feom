#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation

print('Plotting script starts.')

#file_dir = os.environ['PROJECTDIR'] + '/output/'
file_dir = os.getcwd() + '/' # + '/output/'

file_name = sys.argv[1]

calc_name = file_name.split('.')[0]

file_path = file_dir + file_name

med_dir = file_dir + 'media'

pics_dir = file_dir + 'media/pics_' + calc_name

data = np.genfromtxt(file_path)#, skip_header = 2)
qs = data[0,:]
ps = data[1,:]
Ws = data[2:,:]
n_q = len(data)

#**************************************************************
#   Plotting
#**************************************************************
if not os.path.exists(med_dir):
    os.mkdir(med_dir)
    print("Directory " , med_dir ,  " Created ")
else:    
    print("Directory " , med_dir ,  " already exists")

#y_max = 1.3*max(data0[1,1:])
#x_max = 10
#
#if 1==0:
#    print("Plotting images")
#
#    if not os.path.exists(pics_dir):
#        os.mkdir(pics_dir)
#        print("Directory " , pics_dir ,  " Created ")
#    else:    
#        print("Directory " , pics_dir ,  " already exists")
#
#    for t in range(N_t):
#        fig = plt.figure()
#        ax = plt.axes(xlim=(ss[0],ss[-1]), ylim=(0, y_max))
#        #plt.plot(ss, V_s(ss), label="potential")
#        plt.plot(ss, data[t+1,1:], label = "me")
#        plt.plot(ss0, data0[t+1,1:], label = "Stuart")
#        plt.title("t={:5f}".format(ts[t]))
#        #axes = plt.gca()
#        #axes.set_ylim(0,y_max)
#        #axes.set_xlim(-x_max,x_max)
#        plt.legend()
#        #plt.savefig("pics/{:5d}".format(i) + ".png")
#        plt.savefig(pics_dir + "/{:5f}".format(ts[t]) + ".png")
#        plt.clf()
#        plt.close()
#
if 1==1:
    print("Plotting contour")

    #levels = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]
    levels = 10
    fig = plt.figure()
    ax = plt.axes(xlim=(ps[0], ps[-1]), ylim=(qs[0], qs[-1]))
    plt.contour(ps, qs, Ws, levels, cmap='Blues')
    plt.title("Evolution of wavepacket")
    #axes = plt.gca()
    #axes.set_ylim(0,y_max)
    #axes.set_xlim(-x_max,x_max)
    #plt.legend()
    #plt.savefig("pics/{:5d}".format(i) + ".png")
    plt.savefig(med_dir + "/contour_" + calc_name + ".png")
    plt.show()
    plt.clf()
    plt.close()

print("Plotting script finished successfully")
