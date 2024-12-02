#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation

print('Plotting script starts.')

file_name = sys.argv[1]
file_path = file_name
data = np.genfromtxt(file_path)

pics_dir = 'pics'

ts = data[:,0]
N_t = len(ts)
ados = data[:, 1:]

#**************************************************************
#   Plotting
#**************************************************************
y_max = 1.3*max(data0[1,1:])
x_min = 0
x_max = len(ados[0,:])

if 0==1:
    print("Plotting images")

    if not os.path.exists(pics_dir):
        os.mkdir(pics_dir)
        print("Directory " , pics_dir ,  " Created ")
    else:    
        print("Directory " , pics_dir ,  " already exists")

    for t in range(N_t):
        fig = plt.figure()
        ax = plt.axes(xlim=(ss[0],ss[-1]), ylim=(0, y_max))
        #plt.plot(ss, V_s(ss), label="potential")
        plt.plot(ss, data[t+1,1:], label = "me")
        plt.plot(ss0, data0[t+1,1:], label = "Stuart")
        plt.plot(data_tr[t, :], ys, "ro")
        plt.title("t={:5f}".format(ts[t]))
        #axes = plt.gca()
        #axes.set_ylim(0,y_max)
        #axes.set_xlim(-x_max,x_max)
        plt.legend()
        #plt.savefig("pics/{:5d}".format(i) + ".png")
        plt.savefig(pics_dir + "/{:010d}".format(t) + ".png")
        plt.clf()
        plt.close()

if 1==1:
    print("Plotting contour")

    levels = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65]

    fig = plt.figure()
    #ax = plt.axes(xlim=(-x_max, x_max), ylim=(0, y_max))
    plt.contour(ts, ss, np.transpose(data[1:,1:]), levels, cmap='Blues')
    plt.contour(ts, ss0, np.transpose(data0[1:,1:]), levels, cmap='Reds')
    plt.title("Evolution of wavepacket")
    #axes = plt.gca()
    #axes.set_ylim(0,y_max)
    #axes.set_xlim(-x_max,x_max)
    #plt.legend()
    plt.savefig("contour_" + calc_name + ".png")
    #plt.show()
    plt.clf()
    plt.close()

#**************************************************************
#   Animating
#**************************************************************
print(ts[2])
print(ts[1])
print(0.1/abs(ts[2]-ts[1]))
koef=int(round(0.1/abs(ts[2]-ts[1])))
print(koef)
if 1==1 and not tr:
    print("Animating")
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(ss[0], ss[-1]), ylim=(0, y_max))
    #plt.plot(ss, V_s(ss), label="potential")
    line1, = ax.plot([], [], lw=2, label = "data2")
    line2, = ax.plot([], [], lw=2, label = "data1")
    plt.legend()

    # initialization function: plot the background of each frame
    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        return line1, line2,

    # animation function.  This is called sequentially
    def animate(t):
        plt.title("t={:5f}".format(ts[koef*t]))
        print("t={:5f}".format(ts[koef*t]))
        line1.set_data(ss, data[koef*t+1,1:])
        line2.set_data(ss0, data0[koef*t+1,1:])

        return line1, line2,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(N_t/koef), interval=200, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    anim.save('animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
if 1==1 and tr:
    print("Animating")
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(ss[0], ss[-1]), ylim=(0, y_max))
    #plt.plot(ss, V_s(ss), label="potential")
    line1, = ax.plot([], [], lw=2, label = "data2")
    line2, = ax.plot([], [], lw=2, label = "data1")
    trajs, = ax.plot([], [], "ro", label = "trajectories")
    plt.legend()

    # initialization function: plot the background of each frame
    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        trajs.set_data([],[])
        return line1, line2, trajs,

    # animation function.  This is called sequentially
    def animate(t):
        plt.title("t={:5f}".format(ts[koef*t]))
        print("t={:5f}".format(ts[koef*t]))
        line1.set_data(ss, data[koef*t+1,1:])
        line2.set_data(ss0, data0[koef*t+1,1:])
        trajs.set_data(data_tr[koef*t, :], ys)

        return line1, line2, trajs,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(N_t/koef), interval=200, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    #anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    anim.save('animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
print("Plotting script finished successfully")
