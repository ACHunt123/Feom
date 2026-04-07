#!/usr/bin/env python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import json5
import matplotlib
import argparse

#### Parser (for whether to save etc)
parser = argparse.ArgumentParser(description="Plot and optionally save results.")
parser.add_argument('--save', action='store_true', help="save figure as pdf")
parser.add_argument('--pgf', action='store_true', help="use pgf backend for latex. (and save to file)")
args = parser.parse_args()

### Setup the mpl settings 
if args.pgf:
    matplotlib.use("pgf")
plt.rcParams['lines.linewidth'] = 0.6
paper_fontsize = 10
plt.rcParams.update({
    'font.size': paper_fontsize,           # Base font size
    'axes.titlesize': paper_fontsize,      # Title of each subplot
    'axes.labelsize': paper_fontsize+1,      # X/Y axis labels
    'xtick.labelsize': paper_fontsize-4,     # Tick labels
    'ytick.labelsize': paper_fontsize-4,
    'legend.fontsize': paper_fontsize-2 , 
    'figure.titlesize': paper_fontsize + 1 # Bigger if you use suptitle
})
plt.rcParams.update({
    "pgf.texsystem": "pdflatex",  # or xelatex/lualatex
    "font.family": "serif",       # Match LaTeX document
    "text.usetex": True,          # Use LaTeX for all text
    "pgf.rcfonts": False,         # Don't override LaTeX fonts
})

### Which folders to plot (named in the bash script)
folders=["debyeNbead", "debyeLT", "debyePadeN_N","debyeAAA"]
RGplot_folders=["debyeNbead", "debyeLT",'debyePadeN_N','debyeAAA']


### Paths for the files (setup due to the bash script)
path = '{}/beta{}/eta{}/L{}/K{}/dt{}'
endfile = f'{path}/*.out'
FORTfile = f'{path}/tmp/output'
rgfile = f'{path}/*Cothapproximation.txt'

extra=''
### read parameters from batchparams.json
with open("plotparams.json5",'r') as f:
    params = json5.load(f)
Ks = params["Ks"]
Ls = [params["Ls"][-1]]
# Ls = params["Ls"]
dts = params["dts"]
etas = params["etas"]
betas = params["betas"]
beta=betas[0]
gam = params["HRDparams"]["gam"]
Delta = params["HRDparams"]["Delta"]
eps = params["HRDparams"]["eps"]

### Labels for figures
labels={'debyeAAA':"A4",
         'debyePadeN-1_N':"Pade N-1/N", 
        'debyePadeN_N':"Pade [N/N]", 
        'debyeLT':"Mats. + IT",
        'debyeNbead':"Ring-Polymer"
          }

### AXES (vertical Rg, HEOM stacking)
fig, ALLaxes = plt.subplots(2,len(folders), figsize=(6.85,3.5))#, sharex=True, sharey=True)
axes=ALLaxes.flatten()[0:4]
rgaxes=ALLaxes.flatten()[4:8]

### Get all the data (works even if the scripts are still running)
exactfound=False
for folder,ax in zip(folders,axes.flatten()):
    
    for eta in etas:
        for L in Ls:
            for K in Ks:
                for dt in dts:

                    # RG plot [only Pade and AAA write these down. we manually do the others later]
                    if folder in RGplot_folders and L==Ls[0]:
                        rg_file=rgfile.format(folder,beta,eta,L,K,dt)
                        fname = glob.glob(rg_file)[0] if len(glob.glob(rg_file)) > 0 else None
                        if(fname!= None):
                            data = np.loadtxt(fname)
                            w = data[:, 0]
                            signal = data[:, 1]*2/beta
                            rgaxes[RGplot_folders.index(folder)].plot(w,signal,label=f'K={K}')
                            if K==Ks[-1]:#get the exact Rg
                                exactRG = data[:, 2]*2/beta
                        rgaxes[RGplot_folders.index(folder)].set_xlabel(r'$\omega$',labelpad=3)
                        rgaxes[RGplot_folders.index(folder)].set_ylabel(r'$\,\,\,\quad \mathcal{R}^2(\omega)$',labelpad=-1)
                    ################################################################

                    endfile_name =endfile.format(folder,beta,eta,L,K,dt)
                    fname = glob.glob(endfile_name)[0] if len(glob.glob(endfile_name)) > 0 else None
                    if(fname!= None):
                        ### Load the already formatted data from the final output file
                        data = np.loadtxt(fname)
                        t = data[:, 0]
                        signal = data[:, 1]
                    if fname is None:
                        print(f'File {endfile_name} does not exist. looking for the FORTRAN running file')
                        FORTfile_name = FORTfile.format(folder,beta,eta,L,K,dt)
                        fname = glob.glob(FORTfile_name)[0] if len(glob.glob(FORTfile_name)) > 0 else None
                        if fname is None:
                            print(f'File {FORTfile_name} does not exist. Skipping...')
                            continue
                        else:
                            ### calculate from the raw FORTRAN output file
                            data = np.loadtxt(fname)
                            if len(data.shape) == 1:
                                print(f'\nFile {fname} is not formatted correctly. Skipping...\n')
                                continue
                            t= data[:,0]
                            ns=2
                            re_rho = data[:,1:1+ns**2]
                            im_rho = data[:,1+ns**2:]
                            rho = re_rho + 1.j*im_rho
                            rho = rho.reshape((len(t),ns,ns), order='F')  # reshape the data to be a 3D array
                            # format the data to calculate <s_z>, <s_y> and <s_x> and others
                            t = data[:,0]
                            rho11 = rho[:,1,1]  
                            rho00 = rho[:,0,0]
                            rho10 = rho[:,1,0]
                            rho01 = rho[:,0,1]
                            processed_data = np.zeros((len(t),5),dtype=complex)
                            processed_data[:,0] = t
                            processed_data[:,1] = rho11 - rho00  # <s_z>
                            processed_data[:,2] = 1.j*(rho10 - rho01)  # <s_y>
                            processed_data[:,3] = (rho10 + rho01)  # <s_x>   
                            processed_data[:,4] = (1+(rho11 - rho00))/2  # Site 1 population
                            signal= processed_data[:,1]
                        

                    ax.plot(t, signal, label=f'K={K}')

                    ### Get the 'exact'  answer (assumes that AAA converges)
                    if folder=='debyeAAA' and K==Ks[-1] and L==Ls[-1]:
                        # print('\n exact has been found \n')
                        exactfound=True
                        exact_signal=signal
                        exact_t=t
                    #################################################


    ax.set_xlabel(r'$t$',labelpad=1)
    ax.set_ylabel(r'$\left<\sigma_z(t)\right>$',labelpad=3)


    ax.set_title(labels[folder])    

### plot the exact HEOM, and add on the legend
for ax in axes.flatten():
    if exactfound: ax.plot(exact_t,exact_signal, label='Exact',ls='--',zorder=10000,color='k')

### calculate and plot the Rgs for Nbead and LT (we need to calculate these)
for K in Ks:
    hbar=1
    ### Nbead and Nmode
    largenum=50000
    mu = K                  # number of pairs of matsubara modes/ r.p. modes, each pair gives a single exponential term
    N_exp = 1 + mu          # number of exponential terms in the BCF [Temp ind. Exponential, <--- Matsubara Exponentials --->]
    N_mds = 2*mu+1
    betaN = beta/N_mds
    wN=1/(betaN*hbar)
    wns = np.array([2*wN*np.pi*k/N_mds for k in range(0,largenum+1)])
    wks = np.array([2*wN*np.sin(np.pi*k/N_mds) for k in range(0,mu+1)])
        
    ITsum=np.zeros_like(w) #Ishizak-Tanimura sum
    matssum=np.zeros_like(w)
    rpsum=np.zeros_like(w)
    for i in range (1,len(wns)):
        if i<(len(wks)):
            matssum += 1/(w**2+wns[i]**2)
            rpsum += 1/(w**2+wks[i]**2)
            ITsum += 1/(w**2+wns[i]**2)
        else:
            if (wns[i]**2-gam**2)==0:
                print('skipping singularity in IT sum')
            ITsum += 1/(wns[i]**2-gam**2) 
            # ITsum += 1/(wns[i]**2) ### THIS WOULD BE THE 


    mats_modes=matssum*2/beta
    rp_modes=rpsum*2/beta
    matsIT_modes=ITsum*2/beta

    ### plot them (NOTE: hardcoded positions in axes)
    rgaxes[0].plot(w,rp_modes,label=f'K={K}')
    # Rg_ax.plot(w,mats_modes,label='Smooth')
    rgaxes[1].plot(w,matsIT_modes,label=f'K={K}')

### Add on the exact Rg to all plots
for ax in rgaxes.flatten():
    ax.plot(w,exactRG,label=f'Exact',color='k',ls='--')



### Set same y-limits for all axes by Find global min and max y values
all_axes = axes.flatten()  # flatten 2x2 to list
for row,row_axes in enumerate([axes, rgaxes]):
    ymins, ymaxs = [], []
    for ax in row_axes:
        ymin, ymax = ax.get_ylim()
        ymins.append(ymin)
        ymaxs.append(ymax)
    global_ymin = min(ymins)
    global_ymax = max(ymaxs)
    for i,ax in enumerate(row_axes):
        # ax.set_ylim(global_ymin, global_ymax)
        if i!=0:    
            ax.set_yticklabels([])
            ax.set_ylabel('')
        if row==0:
            ax.set_xlim(0, 10)


### add a legend to just the last rgaxes
handles, labels = rgaxes[-1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.tight_layout()
legend=fig.legend(
    by_label.values(),
    by_label.keys(),
    loc='upper center',              # Position at top center of figure
    bbox_to_anchor=(0.9575, 0.6825),      # Fine-tune placement
    # fontsize=7,
    frameon=True,                    # Show frame (border)
    edgecolor='black',              # Set border color
    # title='Legend',                 # Add title
    # title_fontsize=8,               # Optional: title font size
    borderpad=0.1,
    framealpha=None,           # Make the box fully opaque
    facecolor='white'             # Set background color (default is None, i.e. transparent
)
legend.get_frame().set_alpha(1)  # no transparency
for ax in np.concatenate([axes, rgaxes]):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)
legend.get_frame().set_edgecolor('black')
legend.get_frame().set_linewidth(0.5)  # thinner line if needed
legend.get_frame().set_boxstyle('Square')  # force square corners (removes bevel/rounding)

### manually set log scales and limits for the Rg axes
for i,rgax in enumerate(rgaxes):
    rgax.set_xscale('log')
    rgax.set_yscale('log')
    # rgax.set_yscale('symlog', linthresh=10) 
    if i!=0:    
        rgax.set_yticklabels([])
        rgax.set_ylabel('')
    rgax.minorticks_off()
    if beta==500: rgax.set_ylim(1e-1*2/beta,1.1e4*2/beta) #used for beta=500

    rgax.set_xlim(1e-2,100)

### Manual tick labels (HEOM axes)
for i,ax in enumerate(axes):
    ax.set_ylim(-1,1)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticklabels([-1, -0.5, 0, 0.5, 1])
    if i!=0:    
        ax.set_yticklabels([])
        ax.set_ylabel('')

### Padding of whole fig
plt.tight_layout(
    pad=0.,    # space between figure edge and subplots (default: 1.08)
    w_pad=0.05,  # width padding between subplots (default: None -> auto)
    h_pad=0.5   # height padding between subplots (default: None -> auto)
)

### Save everything.
if args.pgf:
    plt.savefig("/home/ach221/software/phd/LaTeX/blob_paper/formatted_paper/figs/results.pgf",bbox_inches='tight')
elif args.save:
    plt.savefig("/home/ach221/software/phd/LaTeX/blob_paper/formatted_paper/figs/results.pdf",bbox_inches='tight')
    plt.show()
else:
    plt.show()
