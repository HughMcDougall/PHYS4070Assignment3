import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

#------------------------------------
filename = "./01"
sparseness      = 1
anim_sparseness = 5
chonk = 2

L=20

mode = 2

#------------------------------------

REAL = np.loadtxt(filename+"-real.dat")[::sparseness]
IMAG = np.loadtxt(filename+"-imag.dat")[::sparseness]
NORM = np.sqrt(REAL**2+IMAG**2)

N = REAL.shape[-1]  #Spatial spacing
M = REAL.shape[0]   #Time spacing

X = np.linspace(-L/2,L/2,N)
T = np.linspace(0,100,M)

#------------------------------------
fig,ax = plt.subplots(1,1, sharex=True, sharey=True)

real, = ax.plot(X, REAL[0,:],lw=1)
imag, = ax.plot(X, IMAG[0,:],lw=1)
norm, = ax.plot(X, NORM[0,:],c='purple',lw=2,)

ax.set_xlim(-L/2,L/2)
ax.set_ylim(-1.5,2.5)
ax.axhline(0,c='k', lw=2, zorder=-1, ls='--')

for i in range(M//anim_sparseness//chonk):
    frameno = i * anim_sparseness * chonk
    ax.plot(X,NORM[frameno,:], alpha=0.03, lw=2 * chonk, c='k')

def animate(i):
    frameno = i * anim_sparseness
    # Update the y-data of the line
    real.set_ydata(REAL[frameno,:])
    imag.set_ydata(IMAG[frameno,:])
    norm.set_ydata(NORM[frameno,:])
    return [real,imag,norm]

ani = animation.FuncAnimation(fig, animate, frames=M//anim_sparseness, interval=1, blit=True)

fig.tight_layout()

if mode == 1:
    #   Phase-ograms
    #------------------------------------

    #Plot phase of mid-point
    fig2, ax2 = plt.subplots(1,1, sharex=True, sharey=True)

    ax2.plot(T,REAL[:,N//2],    lw=1)
    ax2.plot(T,IMAG[:,N//2],    lw=1)
    
    ax2.plot(T,NORM[:,N//2],    c='purple', lw=2)
    ax2.plot(T,-NORM[:,N//2],    c='purple', lw=2)

    ax2.axhline(0,c='k', lw=2, zorder=-1, ls='--')

    #------------------------------------
    fig3, ax3 = plt.subplots(1,1, figsize=(5,5))

    ax3.plot(REAL[:,N//2], IMAG[:,N//2])
    
    ax3.scatter(REAL[:,N//2][0],    IMAG[:,N//2][0],    marker='x')
    ax3.scatter(REAL[:,N//2][-1],   IMAG[:,N//2][-1],   marker='X')
    
    ax3.axis('square')

    #------------------------------------
else:
    #locate indices of peaks
    index_p1 = np.array([np.where(norm == np.max(norm[:N//2]))[0][0] for norm in NORM])
    index_p2 = np.array([np.where(norm == np.max(norm[N//2:]))[0][0] for norm in NORM])

    X_p1, X_p2 = X[index_p1], X[index_p2]
    X_sep = X_p2-X_p1

    maxsep, minsep, avsep  =np.max(X_sep), np.min(X_sep), ((np.max(X_sep) + np.min(X_sep))/2)

    #====
    fig_peak, ax_peak = plt.subplots(1,2,   sharex=False, sharey=True)

    ax_peak[0].plot(X_p1,       T)
    ax_peak[0].plot(X_p2,       T)
    ax_peak[1].plot(X_sep,  T)
    
    ax_peak[0].axvline(0,    c='k', ls='--')
    ax_peak[1].axvline(avsep,    c='k', ls='--')
    ax_peak[1].axvline(minsep,    c='k', ls='--', lw=0.5)
    ax_peak[1].axvline(maxsep,    c='k', ls='--', lw=0.5)

    #====
    
    print("Peak-Peak Separation:")
    print("\t max:\t %.2f" %minsep)
    print("\t min:\t %.2f" %maxsep)
    print("\t  av:\t %.2f" %avsep)
    
plt.show()
