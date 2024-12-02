#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation

print('Plotting script starts.')

pics_dir = 'pics_ado_trend'
file_name = sys.argv[1]
file_path = file_name
data = np.genfromtxt(file_path)
ts = data[:,0]
N_t = len(ts)
xs = list(range(len(data[0,1:])))
ys = data[:,1:]
ys = np.log10(ys)

#**************************************************************
#   Plotting
#**************************************************************
x_min = xs[0]
x_max = xs[-1]
y_min = np.amin(ys)
y_max = np.amax(ys)#*1.1
if m.isnan(y_max) or y_max>50:
    y_max = 30
if m.isnan(y_min) or y_min<50:
    y_min = -10
print("y_min",y_min)
print("y_max",y_max)

if 1==0:
    print("Plotting images")

    if not os.path.exists(pics_dir):
        os.mkdir(pics_dir)
        print("Directory " , pics_dir ,  " Created ")
    else:    
        print("Directory " , pics_dir ,  " already exists")

    for t in range(N_t):
        fig = plt.figure()
        ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min, y_max))
        #plt.yscale('log')
        plt.plot(xs, ys[t,:])
        plt.title("t={:5f}".format(ts[t]))
        plt.ylabel('$\log_{10}$(max magnitude of ADO')
        plt.xlabel('ADO_number')
        #plt.legend()
        plt.savefig(pics_dir + "/{:010d}".format(t) + ".png")
        plt.clf()
        plt.close()

#**************************************************************
#   Animating
#**************************************************************
koef=int(round(0.1/abs(ts[2]-ts[1])))
if 1==1 :
    print("Animating")
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(x_min,x_max), ylim=(y_min, y_max))
    plt.ylabel('log_{10}(max magnitude of ADO')
    plt.xlabel('ADO_number')
    line1, = ax.plot([], [], lw=2, label = "ADO magnitude maxima")
    #plt.legend()

    # initialization function: plot the background of each frame
    def init():
        line1.set_data([], [])
        return line1,

    # animation function.  This is called sequentially
    def animate(t):
        plt.title("t={:5f}".format(ts[koef*t]))
        print("t={:5f}".format(ts[koef*t]))
        line1.set_data(xs, ys[koef*t,:])
        return line1,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(N_t/koef), interval=200, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    anim.save('ado_trend.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
print("Plotting script finished successfully")
