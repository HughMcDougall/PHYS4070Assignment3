'''
This file generates plots for PHYS4070 assignment 3 pt 1
-HM 12/5
'''

#------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import matplotlib.animation as animation

from glob import glob

#------------------------------------
#Load in
FILENAMES = ["./"+a[:-9] for a in glob("*-imag.dat")]

if len(FILENAMES) > 1:
    print("Which output do you want to generate plots for?")
    for i, fname in zip(range(len(FILENAMES)),FILENAMES):
        print("\t%i: \t %s" %(i, fname))
    fi = int(input())
else:
    fi = 0
#------------------------------------
#Load in
filename = FILENAMES[fi]
sparseness      = 1
anim_sparseness = 5
chonk = 2

L    = 20
Tmax = 40

#------------------------------------
# Setup

REAL = np.loadtxt(filename+"-real.dat")[::sparseness]
IMAG = np.loadtxt(filename+"-imag.dat")[::sparseness]
NORM = np.sqrt(REAL**2+IMAG**2)

N = REAL.shape[-1]  #Spatial spacing
M = REAL.shape[0]   #Time spacing

X = np.linspace(-L/2,L/2,N)
T = np.linspace(0,Tmax,M)

# Autodetect mode

#If two peaks in mode 3
if len(find_peaks(NORM[0,:])[0]) == 1:

    #If norm_0 is compact, sech -> mode 2
    var = np.sum(X**2*NORM[0,:]) / np.sum(NORM[0,:])
    if var < 10:
        mode = 2
    else:
        mode = 1
else:
    mode = 3

#------------------------------------
print("Doing mode %i plots for %s" %(mode,filename))

#------------------------------------
# ANIMATION
fig,ax = plt.subplots(1,1, sharex=True, sharey=True, figsize=(6,3))

real, = ax.plot(X, REAL[0,:],   lw=1,   alpha=0.5, label = 'Re ($\psi$)')
imag, = ax.plot(X, IMAG[0,:],   lw=1,   alpha=0.5, label = 'Im ($\psi$)')
norm, = ax.plot(X, NORM[0,:],c='purple',lw=2, label = "$|(\psi)|$")

ax.set_xlim(-L/2,L/2)
ax.set_ylim(-1.5,2.5)
ax.axhline(0,c='k', lw=2, zorder=-1, ls='--')
ax.legend(loc='best')

for i in range(M//anim_sparseness//chonk):
    frameno = i * anim_sparseness * chonk
    ax.plot(X,NORM[frameno,:], alpha=0.01, lw=2 * chonk, c='k')

def animate(i):
    frameno = i * anim_sparseness
    # Update the y-data of the line
    real.set_ydata(REAL[frameno,:])
    imag.set_ydata(IMAG[frameno,:])
    norm.set_ydata(NORM[frameno,:])
    return [real,imag,norm]

ani = animation.FuncAnimation(fig, animate, frames=M//anim_sparseness, interval=1, blit=True)

fig.tight_layout()

if mode == 1 or mode==2:
    #   Phase-ograms
    #------------------------------------

    if mode == 1:
        index_p = np.array([       N//2 for norm in NORM])
    elif mode == 2:
        index_p = np.array([       np.where(norm == np.max(norm))[0][0] for norm in NORM])
    X_p   = X[index_p]
    X_av  = np.array([np.sum(X*norm) / np.sum(norm) for norm in NORM])

    #------------------------------------

    #Plot phase of peak or midpoint
    fig2, ax2 = plt.subplots(1,1, sharex=True, sharey=True, figsize = (6,3))

    realpeak = np.array([real[i] for real, i in zip(REAL,index_p)])
    imagpeak = np.array([imag[i] for imag, i in zip(IMAG,index_p)])
    normpeak = np.sqrt(realpeak**2+imagpeak**2)

    real_av = np.array([np.sum(real * norm) / np.sum(norm) for real,norm in zip(REAL,NORM)])
    imag_av = np.array([np.sum(imag * norm) / np.sum(norm) for imag,norm in zip(IMAG,NORM)])
    norm_av = np.array([np.sum(norm * norm) / np.sum(norm) for norm,norm in zip(NORM,NORM)])

    if mode==2:
        #Plot Position of peak
        fig_peakloc, ax_peakloc = plt.subplots(1,1, sharex=True, sharey=True, figsize = (6,2))

        ax_peakloc.plot(T, X_p, c= 'dodgerblue', lw=2, label="Peak Location")
        ax_peakloc.plot(T, X_av, c = 'blue', lw=1, label="Average Position")
        ax_peakloc.legend()
        fig_peakloc.supxlabel("Time")
        ax_peakloc.axhline(0, c='k', lw=2, zorder=-1, ls='--')
        fig_peakloc.tight_layout()

    #----------------------------
    ax2.plot(T,realpeak,    lw=1)
    ax2.plot(T,imagpeak,    lw=1)
    
    ax2.plot(T, normpeak,    c='purple', lw=2)
    ax2.plot(T,-normpeak,    c='purple', lw=2)

    ax2.axhline(0, c='k', lw=2, zorder=-1, ls='--')

    #------------------------------------
    # Phase diagram
    fig3, ax3 = plt.subplots(1,1, figsize=(3,3))

    thet = np.arctan2(imagpeak, realpeak)
    thet = np.linspace(np.min(thet),np.max(thet),256)

    ax3.plot(realpeak, imagpeak)
    ax3.plot(np.cos(thet)*np.max(normpeak), np.sin(thet)*np.max(normpeak),c='k',lw=1)
    ax3.plot(np.cos(thet)*np.min(normpeak), np.sin(thet)*np.min(normpeak),c='k',lw=1)

    ax3.plot(real_av, imag_av)
    
    ax3.set_xlabel('Re ($\psi$)')
    ax3.set_ylabel('Im ($\psi$)')
    
    ax3.axis('square')
    fig3.tight_layout()

    #------------------------------------
elif mode==3:
    #locate indices of peaks
    index_p1 = np.array([       np.where(norm[:N//2] == np.max(norm[:N//2]))[0][0] for norm in NORM])
    index_p2 = np.array([N//2 + np.where(norm[N//2:] == np.max(norm[N//2:]))[0][0] for norm in NORM])

    X_p1, X_p2 = X[index_p1], X[index_p2]
    X_sep      = X_p2 - X_p1

    X_p1_av  = np.array([np.sum(X[:N//2]*norm[:N//2]) / np.sum(norm[:N//2]) for norm in NORM])
    X_p2_av  = np.array([np.sum(X[N//2:]*norm[N//2:]) / np.sum(norm[N//2:]) for norm in NORM])
    X_sep_av = X_p2_av - X_p1_av


    maxsep, minsep, avsep  =np.max(X_sep), np.min(X_sep), ((np.max(X_sep) + np.min(X_sep))/2)

    #====
    # Peak-Peak separations
    fig_peak, ax_peak = plt.subplots(1,2,   sharex=False, sharey=True, figsize=(6,3))

    ax_peak[0].plot(X_p1,       T)
    ax_peak[0].plot(X_p2,       T)
    ax_peak[1].plot(X_sep,  T)

    ax_peak[0].plot(X_p1_av,       T, lw=1)
    ax_peak[0].plot(X_p2_av,       T, lw=1)
    ax_peak[1].plot(X_sep_av,      T, lw=1)
    
    ax_peak[0].axvline(0,    c='k', ls='-')
    
    ax_peak[1].axvline(avsep,    c='k', ls='--')
    ax_peak[1].axvline(minsep,    c='k', ls='--', lw=0.5)
    ax_peak[1].axvline(maxsep,    c='k', ls='--', lw=0.5)

    #ax_peak[0].grid()
    #ax_peak[1].grid()


    ax_peak[0].set_xlim(-L/2,L/2)
    ax_peak[1].set_xlim(0,L)
    ax_peak[0].set_ylim(0,Tmax)

    ax_peak[0].set_title("Peak Locations")
    ax_peak[1].set_title("Peak Separation")
    
    fig_peak.supylabel("Time")
    fig_peak.supxlabel("Peak Separation")
    fig_peak.tight_layout()

    #====
    
    print("Peak-Peak Separation:")
    print("\t max:\t %.2f" %maxsep)
    print("\t min:\t %.2f" %minsep)
    print("\t  av:\t %.2f" %avsep)

#==========================
plt.figure(figsize = (8,4))
dx = L/(N-1)
NORMALL = np.sum(NORM,axis=1)* dx

plt.plot(T,NORMALL)
plt.xlabel("Time")
plt.ylabel("Total Signal Norm, $\int |\psi|^2 dx$")
plt.ylim(0,np.max(NORMALL*1.2))
plt.show()
