import numpy as np
import matplotlib.pyplot as plt
from glob import glob
#------------------------------
'''
The following two lines should narrow the files for plotting down to only those that
currently exist in the ./results folder. If you're having issues with this, remove line
11 and fill the TARGS vector with the files manually.
'''
TARGS = ['energies-4.dat', 'energies-5.dat', 'energies-6.dat', 'energies-8.dat', 'energies-10.dat']
TARGS = [targ[2:] for targ in glob("./*") if targ[2:] in TARGS]
TARGS = []

if len(TARGS)==0:
    print("There seems to be an error loading files for plotting. please open the _Q2plot.py file for help")
    input("Press enter to exit")
    assert(False)
#------------------------------
def dE(x,thet):
    return np.sqrt(1-x*np.sin(thet)**2)

def E(x, N_int=128):
    dthet = np.pi / 2 / (N_int+1)
    out = 0

    for i in range(N_int):
        out+= dthet * dE(x, dthet * (i+1/2))
    return(out)

def eps(g, N_int=128):
    out = -1/2 / np.pi * (2+g) * E( 8 * g / (2+g)**2 , N_int)
    return(out)

def ddeps(g,dg = 0.0001, N_int=128):
    if g==0: return( ddeps(g+dg, N_int = N_int) )

    out = (eps(g+dg, N_int) - 2*eps(g, N_int) + eps(g-dg, N_int)) / dg**2
    return(out)

#------------------------------
# Analytical Approximations
Nplot = 128
dgplot = 0.0001

Gplot   = np.linspace(0,8,Nplot)
Enplot   = np.array([eps(g, N_int = 128) for g in Gplot])
DDEnplot = np.array([ddeps(g,dgplot, N_int = 512) for g in Gplot])

#------------------------------
# Make plots, plot Analytical
fig_en, ax_en = plt.subplots(1,1, figsize=(6,3))
fig_dden, ax_dden = plt.subplots(1,1,figsize=(6,3))

ax_en.plot(Gplot,Enplot, c = 'k', label = 'Analytical $N = \infty$ Solution')
ax_dden.plot(Gplot,DDEnplot, c = 'k', label = 'Analytical $N = \infty$ Solution')
#------------------------------
# Plot Simulated
for targ in TARGS:

    #Get targ
    name = "N="+targ[9:-4]
    X = np.loadtxt(targ)

    #===
    #Extract data
    G       = X[:,0]
    En      = X[:,1]
    En_norm = X[:,2]

    N = len(G)

    #===
    #Calculate derivatives
    ddEn_norm = np.zeros(N)
    dg = (G[-1] - G[0]) / (N-1)
    for i in range(1,N-1):
        ddEn_norm[i] = (En_norm[i+1] - 2 * En_norm[i] + En_norm[i-1]) / dg **2
    ddEn_norm[0]=ddEn_norm[1]
    ddEn_norm[-1]=ddEn_norm[-2]

    #===
    #Do plots
    ax_en.plot(G,En_norm,marker = 'x', lw=0)
    ax_dden.plot(G,ddEn_norm,marker = 'x', label=name)

#------------------------------
# Labels etc
ax_en.axhline(-0.5, ls='--', c='k', lw=1)
ax_en.axhline(0.0, ls='-', c='k', lw=1)
ax_en.axline([0,0],slope=-0.25, ls='--', c='k', lw=1)
ax_dden.axhline(0.0, ls='-', c='k', lw=1)

ax_en.set_xlim(0,8)
ax_dden.set_xlim(0,8)

ax_en.legend()
ax_dden.legend()

fig_en.supxlabel("Strength of Spin-Spin Interaction, g")
fig_dden.supxlabel("Strength of Spin-Spin Interaction, g")

fig_en.supylabel("$\epsilon_0$ / N")
fig_dden.supylabel("$\partial^2 \epsilon_0 / \partial g^2 $ / N")

fig_en.tight_layout()
fig_dden.tight_layout()

plt.show()
