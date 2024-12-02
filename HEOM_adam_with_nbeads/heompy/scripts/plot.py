#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation

contour = True
averages = True
pics = True
animate = True


#contour = False
#averages = False
#pics = False
#animate = False

p_hist = True
p_tr = True

#p_hist = False
#p_tr = False

print('Plotting script starts.')

file_list = []
files_trajs = []
dens = True
for ind in range(len(sys.argv)-1):
    if dens and sys.argv[ind+1]!="x":
        file_list.append(sys.argv[ind+1])
    if sys.argv[ind+1]=="x":
        dens = False
    if not dens and sys.argv[ind+1]!="x":
        files_trajs.append(sys.argv[ind+1])

print("Plotting density files:")
for item in file_list:
    print(item)
print("Plotting trajectory files:")
for item in files_trajs:
    print(item)
print("Loading and processing data")
tr=True
if len(files_trajs) == 0:
    tr = False

pics_dir = 'pics'

# Reading in data
data_list_raw = []
for item in file_list:
    data_list_raw.append(np.genfromtxt(item))

qs = data_list_raw[0][0,1:]
ts = data_list_raw[0][1:,0]
N_t = len(ts)
N_f = len(data_list_raw)
N_trf = len(files_trajs)
# Adjusting all data to first file's time grid
data_list = []
qs_list = []
avg_list = []
data_list.append(data_list_raw[0])
qs_list.append(qs)

avg0 = np.zeros([N_t,2])
for t in range(N_t):
    avg0[t,0] = data_list_raw[0][t+1,0]
    avg0[t,1] = np.average(qs, weights=data_list_raw[0][t+1,1:])
avg_list.append(avg0)

for datafile in data_list_raw[1:]:
    N_s0 = datafile.shape[1]-1
    qs0 = datafile[0, 1:]
    ts0 = datafile[1:,0]
    N_t0 = len(ts0)

    data0 = np.zeros([N_t+1,N_s0+1])
    data0[0,0] = 0.0
    data0[0,1:] = qs0
    avg0 = np.zeros([N_t,2])
    for t in range(N_t):
        tx = np.searchsorted(ts0,ts[t])
        if tx == N_t0:
            tx = N_t0-1
        data0[t+1,:] = datafile[tx+1,:]
        avg0[t,1] = np.average(qs0, weights=data0[t+1,1:])
        avg0[t,0] = data0[t+1,0]
    data_list.append(data0)
    qs_list.append(qs0)
    avg_list.append(avg0)

# Trajectories
if tr:
    N_bins = 20
    coef_norm = 0.9
    traj_list = []
    ys_list = []
    hist_list = []
    tr_avg_list = []

    for f in range(N_trf):
        data_trajs = np.genfromtxt(files_trajs[f])
        N_trajs = data_trajs.shape[1]-1
        t_trajs = data_trajs[0:,0]
        N_t_trajs = len(t_trajs)
        data_tr = np.zeros([N_t, N_trajs+1])
        avg0 = np.zeros([N_t, 2])
        ys = np.linspace(0,0.2,N_trajs)
        for t in range(N_t):
            tx = np.searchsorted(t_trajs,ts[t])
            if tx == N_t_trajs:
                tx = N_t_trajs-1
            data_tr[t,:] = data_trajs[tx,:]
        # Generating histogram data
        data_hist =np.zeros((N_t,N_bins+1,2))
        for t in range(N_t):
            data_hist[t,1:,1], bin_edges = np.histogram(data_tr[t, 1:], bins=N_bins, density = True)
            data_hist[t,0,1] = data_tr[t,0]
            for j in range(N_bins):
                data_hist[t,j+1,0] = (bin_edges[j+1]+bin_edges[j])/2
        data_hist[:,1:,1] = coef_norm*data_hist[:,1:,1]
        avg0[:,0] = data_tr[:,0]
        avg0[:,1] = np.average(data_tr[:, 1:], axis=1)

        traj_list.append(data_tr)
        hist_list.append(data_hist)
        tr_avg_list.append(avg0)
        ys_list.append(ys)


#**************************************************************
#   Plotting
#**************************************************************
x_min = qs[0]
x_max = qs[-1]
y_max = 1.3*np.amax(data_list[0][1,1:])
y_min = 0.0
maxval = np.amax(avg_list[0][:,1])
minval = np.amin(avg_list[0][:,1])
if maxval>15.0 or m.isnan(maxval):
    maxval = 4.0
if minval<-15.0 or m.isnan(minval):
    minval = -1.0
rang = maxval - minval
avg_max = maxval + 0.1 * rang
avg_min = minval - 0.1 * rang
koef=int(round(0.1/abs(ts[2]-ts[1])))

#******************************************************************************
cols = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
cols = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
cnt = 0
cmaps = ['Blues', 'Reds', 'Greens', 'Greys', 'Purples', 'Oranges', 'YlOrBr',
        'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu',
        'PuBuGn', 'BuGn', 'YlGn']
if contour:
    print("Plotting contour")

    levels = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65]

    fig = plt.figure()
    for f in range(N_f):
        plt.contour(ts, qs_list[f], np.transpose(data_list[f][1:,1:]),
                    levels, cmap=cmaps[f])
    # Legend for contour plots is difficult
    plt.title("Evolution of wavepacket")
    plt.savefig("contour.png")
    #plt.show()
    plt.clf()
    plt.close()
#******************************************************************************
if averages:
    lw = 1.0
    fig = plt.figure(figsize=(8, 6), dpi=160)
    ax = plt.axes(xlim=(ts[0],ts[-1]), ylim=(avg_min, avg_max))
    # Plotting densities
    for f in range(N_f):
        plt.plot(avg_list[f][:,0], avg_list[f][:,1],
                                        label=file_list[f], c=cols[cnt], lw=lw)
        cnt += 1
    # Plotting trajectories
    if tr:
        for f in range(N_trf):
            plt.plot(tr_avg_list[f][:,0], tr_avg_list[f][:,1],
                                    label=files_trajs[f], c=cols[cnt], lw=lw)
            cnt+=1
    cnt = 0
    plt.title("Evolution of averages over time")
    plt.xlabel("t/a.u.")
    plt.ylabel("q/a.u.")
    plt.legend(loc=0)
    plt.savefig("averages.png")

    #plt.show()

    plt.clf()
    plt.close()

#******************************************************************************

if pics:
    print("Plotting images")
    # Creating folder
    if not os.path.exists(pics_dir):
        os.mkdir(pics_dir)
        print("Directory " , pics_dir ,  " Created ")
    else:    
        print("Directory " , pics_dir ,  " already exists")
    # Plotting
    for t in range(0,int(N_t), koef):
        print("Image ", int(t/koef), " out of ", int(N_t/koef))
        fig = plt.figure()
        ax = plt.axes(xlim=(x_min, x_max), ylim=(y_min, y_max))
        # Plotting densities
        for f in range(N_f):
            plt.plot(qs_list[f], data_list[f][t+1,1:], label = file_list[f],
                                                                c=cols[cnt])
            plt.axvline(avg_list[f][t,1], c=cols[cnt])
            cnt += 1
        # Plotting trajectories
        if tr:
            for f in range(N_trf):
                if p_hist:
                    plt.plot(hist_list[f][t,1:,0],hist_list[f][t,1:,1],
                                                label=files_trajs[f], c=cols[cnt])
                if p_tr:
                    plt.plot(traj_list[f][t, 1:], ys_list[f], "o", markersize=1,
                                                            alpha=0.1, c=cols[cnt])
                plt.axvline(tr_avg_list[f][t,1], c=cols[cnt])
                cnt+=1
        plt.title("t={:3f}".format(ts[t]))
        plt.legend(loc=1)
        plt.savefig(pics_dir + "/{:05d}".format(int(t/koef)) + ".png")
        plt.clf()
        plt.close()
        cnt = 0

#**************************************************************
#   Animating
#**************************************************************
koef=int(round(0.1/abs(ts[2]-ts[1])))
#******************************************************************************
if animate:
    print("Animating")
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(x_min, x_max), ylim=(y_min, y_max))
    lines = []
    lines_avg = []
    lines_hist = []
    lines_trajs = []
    lines_tr_avg = []
    for f in range(N_f):
        line = ax.plot([], [], lw=2, label = file_list[f], c=cols[cnt])[0]
        avg = ax.axvline(x=0.0, c=cols[cnt])
        cnt +=1
        lines.append(line)
        lines_avg.append(avg)
    if tr:
        for f in range(N_trf):
            trajs = ax.plot([], [], "o", markersize=1, alpha=0.1, c=cols[cnt])[0]
            hist = ax.plot([], [], lw=2, label=files_trajs[f], c=cols[cnt])[0]
            tr_avg = ax.axvline(x=0.0, c=cols[cnt])
            cnt +=1
            lines_trajs.append(trajs)
            lines_hist.append(hist)
            lines_tr_avg.append(tr_avg)

    cnt = 0
    plt.legend(loc=1)

    # initialization function: plot the background of each frame
    def init():
        for line in lines:
            line.set_data([], [])
        if tr:
            for line in lines_hist:
                line.set_data([], [])
            for line in lines_trajs:
                line.set_data([], [])
        return lines

# animation function.  This is called sequentially
    def animate(t):
        plt.title("t={:2f}".format(ts[koef*t]))
        print("t={:5f}".format(ts[koef*t]))
        for f in range(N_f):
            lines[f].set_data(qs_list[f], data_list[f][koef*t+1,1:])
            lines_avg[f].set_xdata(avg_list[f][koef*t,1])
        if tr:
            for f in range(N_trf):
                if p_hist:
                    lines_hist[f].set_data(hist_list[f][koef*t,1:,0],
                                                        data_hist[koef*t,1:,1])
                if p_tr:
                    lines_trajs[f].set_data(traj_list[f][koef*t, 1:], ys_list[f])
                lines_tr_avg[f].set_xdata(tr_avg_list[f][koef*t,1])
        return lines

# call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(N_t/koef), interval=200, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
    anim.save("animation.mp4", fps=30, extra_args=['-vcodec', 'libx264'])

print("Plotting script finished successfully")
