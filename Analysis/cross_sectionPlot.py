import mytools as mt
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

from scipy.interpolate import spline

font = {'family':'DejaVu Sans',
        'size': 20}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

G1n = mt.ENDF_crossection("/Users/jeffreyburggraf/Downloads/1n.txt")
G2n = mt.ENDF_crossection("/Users/jeffreyburggraf/Downloads/2n.txt")
GFiss = mt.ENDF_crossection("/Users/jeffreyburggraf/Downloads/fiss.txt")


def get_m(x):

    result = np.where(np.array(x)<10.5, 1.61*np.e**(-0.54*np.array(x)),0)

    return result

def r_one(x):
    return np.ones_like(x)

fig, axs = plt.subplots(2,1,sharex=True, figsize=(10,6))
axs = axs.flatten()

plt.subplots_adjust(bottom=0.13)


for i,m in enumerate([r_one, get_m]):
    ax = axs[i]
    ax.plot(G1n.x, G1n.y*m(G1n.x),linestyle="--", label="$(\gamma, n)$", color="black")
    if i == 0:
        ax.plot(G2n.x, G2n.y*m(G2n.x), linestyle="-.",label="$(\gamma, 2n)$", color="black")

    ax.plot(GFiss.x, GFiss.y*m(GFiss.x), label="$(\gamma, fiss)$", color="black")
    if i == 0:
        ax.legend()
        ax.set_ylim(0, 1.15*max(G1n.y))

    ax.set_xticks(np.arange(5,15,1))
    ax.grid()
    ax.set_xlim(5, 15)
    ax.ticklabel_format(style='sci', axis="y", scilimits=(0, 0))

    ax.set_ylabel("Cross Section (b)" if i == 0 else r"(Cross Section)$\times p(E)$")
    if i==1:
        ax.set_xlabel("Incident photon energy [MeV]")


        print("Fiss: {0:.2E}, 1n:{1:.2E}".format(sum( GFiss.y*m(GFiss.x)),sum( G1n.y*m(G1n.x))))

plt.show()