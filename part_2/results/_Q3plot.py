import numpy as np
import matplotlib.pyplot as plt
#------------------------------
# Get Data
targname = "./quenchsim"

REAL = np.loadtxt(targname + "-real.dat")[:-1]
IMAG = np.loadtxt(targname + "-imag.dat")[:-1]
NORM = np.sqrt(REAL**2 + IMAG**2)

SZ  = np.loadtxt(targname + "-Sz.dat")
SX  = np.loadtxt(targname + "-Sx.dat")
CXX = np.loadtxt(targname + "-Cxx.dat")

M = NORM.shape[0]
L = NORM.shape[1]
N = int(np.log2(L))

T = np.loadtxt(targname + "-T.dat")

#------------------------------
# Plot 1- System Initial State

fig_initstate, ax_initstate = plt.subplots(2,1, sharex=True, sharey=True, figsize=(6,3))

    
for i in range(M):
    norm = np.insert(NORM[i,:],0,0)
    real = np.insert(REAL[i,:],0,0)
    imag = np.insert(IMAG[i,:],0,0)

    if i==0:
        ax_initstate[0].step(x=np.linspace(-1,L,L+1), y=norm, c='k', lw=1, alpha=1.00, zorder=100)
        ax_initstate[1].step(x=np.linspace(-1,L,L+1), y=real, c='k', lw=1, alpha=1.00, zorder=100)
        ax_initstate[1].step(x=np.linspace(-1,L,L+1), y=imag, c = 'k', lw=1, alpha=1.0, zorder=100)
    else:    
        ax_initstate[0].step(x=np.linspace(-1,L,L+1), y=norm, c='purple', lw=0.1, alpha=0.5)
        ax_initstate[1].step(x=np.linspace(-1,L,L+1), y=real, c='dodgerblue', lw=0.1, alpha=0.25)
        ax_initstate[1].step(x=np.linspace(-1,L,L+1), y=imag, c = 'orange', lw=0.1, alpha=0.25)
    
ax_initstate[0].set_xlim(-1,50)
ax_initstate[0].set_ylabel("Norm", fontsize = 12)
ax_initstate[1].set_ylabel("Real / Imag", fontsize = 12)

fig_initstate.supxlabel("State Number")
fig_initstate.tight_layout()

#------------------------------
# Plot 2- evolution of system state

Nstates = 6
fig, ax = plt.subplots(2,Nstates, sharex=True, sharey=True, figsize=(6,3))
for i in range(Nstates):
    ax[0][i].plot(T,NORM[:,i], c='purple' )
    ax[1][i].plot(T,REAL[:,i], c='dodgerblue' )
    ax[1][i].plot(T,IMAG[:,i], c='orange' )

    ax[0][i].axhline(0, ls = '--', c='k', zorder=-10)
    ax[1][i].axhline(0, ls = '--', c='k', zorder=-10)

    ax[0][i].set_title("State %i" %i)

ax[0][0].set_ylabel("Norm", fontsize = 12)
ax[1][0].set_ylabel("Real / Imag", fontsize = 12)


fig.supxlabel("Time")
fig.tight_layout()

#------------------------------
# Plot 3- evolution of Observables

fig_spins, ax_spins = plt.subplots(3,1, sharex=True, sharey=True, figsize=(6,3))

ax_spins[0].plot(T,SZ/N, c='crimson')
ax_spins[1].plot(T,SX/N, c='gold')
ax_spins[2].plot(T,CXX/N, c='navy')

ax_spins[0].set_ylabel("$S_z / N$", fontsize = 12)
ax_spins[1].set_ylabel("$S_x / N$", fontsize = 12)
ax_spins[2].set_ylabel("$C_{xx} / N$", fontsize = 12)

ax_spins[0].axhline(0, ls = '--', c='k', zorder=-10)
ax_spins[1].axhline(0, ls = '--', c='k', zorder=-10)
ax_spins[2].axhline(0, ls = '--', c='k', zorder=-10)

fig_spins.supxlabel("Time")
fig_spins.tight_layout()
#------------------------------


plt.show()
