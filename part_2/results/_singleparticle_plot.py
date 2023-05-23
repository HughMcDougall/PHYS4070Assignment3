import numpy as np
import matplotlib.pyplot as plt
#------------------------------
targname = "./quenchsim"
targname = "./quench_minor"
sparse = 5
extent = int(-1)

# Get Data
REAL = np.loadtxt(targname + "-real.dat")[:extent:sparse]
IMAG = np.loadtxt(targname + "-imag.dat")[:extent:sparse]
NORM = np.sqrt(REAL**2 + IMAG**2)

SZ  = np.loadtxt(targname + "-Sz.dat")[:extent:sparse]
SX  = np.loadtxt(targname + "-Sx.dat")[:extent:sparse]
CXX = np.loadtxt(targname + "-Cxx.dat")[:extent:sparse]

M = NORM.shape[0]
L = NORM.shape[1]
N = int(np.log2(L))

T = np.loadtxt(targname + "-T.dat")[:extent:sparse]

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
    ax[0][i].plot(NORM[:,i], c='purple' )
    ax[1][i].plot(REAL[:,i], c='dodgerblue' )
    ax[1][i].plot(IMAG[:,i], c='orange' )

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

ax_spins[0].set_ylabel("$S_z$", fontsize = 12)
ax_spins[1].set_ylabel("$S_x$", fontsize = 12)
ax_spins[2].set_ylabel("$C_{xx}$", fontsize = 12)

ax_spins[0].axhline(0, ls = '--', c='k', zorder=-10)
ax_spins[1].axhline(0, ls = '--', c='k', zorder=-10)
ax_spins[2].axhline(0, ls = '--', c='k', zorder=-10)

fig_spins.supxlabel("Time")
fig_spins.tight_layout()
#------------------------------

# Plot 4- Collapsing to individual particle
def int_to_binarray(x, m=None):
    
    X = bin(x)[2:]
    N = len(X)
    
    if type(m)==type(None): m=N
    assert m>=N, "bad length"
    
    out=np.zeros(m, dtype='int32')
    for i in range(N):
        out[-i-1] = int(X[-i-1])
    return out

A = np.stack([int_to_binarray(x, N) for x in range(L)],axis=1)
A = A[:,::-1]   #Prevalence of Z spin up
B = 1-A         #Prevalence of Z spin down

#------------------------------
# z- spin states
REAL_PART_POS = np.dot(A,REAL.T).T[:,0]
IMAG_PART_POS = np.dot(A,IMAG.T).T[:,0]
NORM_PART_POS = np.dot(A,NORM.T).T[:,0]
#NORM_PART_POS = np.sqrt(REAL_PART_POS**2 + IMAG_PART_POS**2)

REAL_PART_NEG = np.dot(B,REAL.T).T[:,0]
IMAG_PART_NEG = np.dot(B,IMAG.T).T[:,0]
NORM_PART_NEG = np.dot(B,NORM.T).T[:,0]
#NORM_PART_NEG = np.sqrt(REAL_PART_NEG**2 + IMAG_PART_NEG**2)

NORM_PART = np.sqrt(NORM_PART_POS**2 + NORM_PART_NEG**2)
PHASE_1 = np.arctan2(IMAG_PART_POS , REAL_PART_POS)
PHASE_2 = np.arctan2(IMAG_PART_NEG , REAL_PART_NEG)
PHASE = PHASE_2 - PHASE_1

#------------------------------
# x - spin states
REAL_PART_X_POS = (REAL_PART_POS + REAL_PART_NEG) / np.sqrt(2)
IMAG_PART_X_POS = (IMAG_PART_POS + IMAG_PART_NEG) / np.sqrt(2)
REAL_PART_X_NEG = (REAL_PART_POS - REAL_PART_NEG) / np.sqrt(2)
IMAG_PART_X_NEG = (IMAG_PART_POS - IMAG_PART_NEG) / np.sqrt(2)

PHASE_1_X = np.arctan2(IMAG_PART_X_POS , REAL_PART_X_POS)
PHASE_2_X = np.arctan2(IMAG_PART_X_NEG , REAL_PART_X_NEG)

NORM_PART_X_POS = np.sqrt(REAL_PART_X_POS**2 + IMAG_PART_X_POS**2)
NORM_PART_X_NEG = np.sqrt(REAL_PART_X_NEG**2 + IMAG_PART_X_NEG**2)
#------------------------------
#Net Spins

SPIN_Z = NORM_PART_POS - NORM_PART_NEG
SPIN_X = NORM_PART_X_POS - NORM_PART_X_NEG
#------------------------------

fig_bypart, ax_bypart = plt.subplots(1,2,sharex=True, sharey=True, figsize=(6,3.5))
ax_bypart[0].plot(NORM_PART_POS, NORM_PART_NEG,lw=1/sparse)
ax_bypart[1].plot(NORM_PART_X_POS, NORM_PART_X_NEG,lw=1/sparse)

ax_bypart[0].set_box_aspect(1)
ax_bypart[1].set_box_aspect(1)


ax_bypart[0].set_title('Z-Spins')
ax_bypart[1].set_title('X-Spins')
ax_bypart[0].set_xlim(-2,2)
ax_bypart[0].set_ylim(-2,2)

fig_bypart.supxlabel("Up")
fig_bypart.supylabel("Down")

#------------------------------
fig_bypart_xz, ax_bypart_xz = plt.subplots(1,2,sharex=True, sharey=True, figsize=(6,3.5))
ax_bypart_xz[0].plot(REAL_PART_POS, IMAG_PART_POS,lw=1/sparse)
ax_bypart_xz[0].plot(REAL_PART_NEG, IMAG_PART_NEG,lw=1/sparse)
ax_bypart_xz[1].plot(REAL_PART_X_POS, IMAG_PART_X_POS,lw=1/sparse)
ax_bypart_xz[1].plot(REAL_PART_X_NEG, IMAG_PART_X_NEG,lw=1/sparse)

ax_bypart_xz[0].set_aspect('equal')
ax_bypart_xz[1].set_aspect('equal')

ax_bypart_xz[0].set_title('Z-Spin Modes')
ax_bypart_xz[1].set_title('X-Spin Modes')

fig_bypart_xz.supxlabel("Re")
fig_bypart_xz.supylabel("Im")
#------------------------------
fig_spindir, ax_spindir = plt.subplots(1,1, figsize=(4,4))
ax_spindir.plot(SPIN_X, SPIN_Z,lw=1)
ax_spindir.set_aspect('equal')
ax_spindir.set_box_aspect(1)
fig_spindir.supxlabel("Net X-Spin")
fig_spindir.supylabel("Net Z-Spin")
fig_spindir.tight_layout()
#------------------------------

plt.show()
