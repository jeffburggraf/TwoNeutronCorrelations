import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
font = {'family': 'DejaVu Sans',
        'size': 28}
mpl.rc('font', **font)
mpl.rc("savefig", dpi=400)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command

mpl.rc('text', usetex=True)


np.random.seed(2)
# centers = [26, 47, 73, 96, 120]
centers = np.arange(24,24*6,24)
amps = [34, 25, 15, 7.5, 3]
sigs = [13, 7,  7,  5 , 8]
cross_offsets = [0,0,0,0,0]
cross_amp_offsets = [0.65,0.8,2,4,5.5]


def vectorAngle(v1,v2,radians=True,eps=0.00001):
    dot=sum([x1*x2 for x1,x2 in zip(v1,v2)])
    mag1mag2=np.sqrt(sum([x*x for x in v1]))*np.sqrt(sum([x*x for x in v2]))

    cos=dot/(mag1mag2)
    if cos>1.-eps or cos<-1.+eps:
        return 0
    else:
        return np.arccos(cos) if radians else (180/np.pi)*np.arccos(cos)

data_corr = []
data_cross = []

for c, a, sig, cross_offset, cross_amp_in in zip(centers, amps, sigs, cross_offsets, cross_amp_offsets):
    th = np.pi/180.*c
    x = 1.25*np.cos(th)
    y = 1.25*np.sin(th)
    # plt.plot([c]*100,np.linspace(0,14000, 100))
    angles_c = []
    angles_cross = []
    l = int(a*100*sig)
    for i in range(int(l)):
        v1_c = (1.25,0,np.random.rand()*0.381)
        v2_c = (x,y,np.random.rand()*0.381)
        z_cross = np.random.rand() * 0.381
        v1_cross = (1.25, 0, z_cross)
        v2_cross = (x, y,z_cross + 0.1*np.random.randn())
        angle_c = vectorAngle(v1_c, v2_c,radians=False)
        if i<cross_amp_in*l*0.04:
            angle_cross = vectorAngle(v1_cross, v2_cross,radians=False) + cross_offset
            angles_cross.append(angle_cross)

        angles_c.append(angle_c)

    data_cross.extend(np.array(angles_cross) + np.random.randn(len(angles_cross))*sig*0.2)
    data_corr.extend(np.array(angles_c) + np.random.randn(l)*sig*0.2)

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

########################
binwidth = 2
#####################
plt.hist(data_corr, np.arange(20,135,binwidth), histtype='step', linewidth=3, label='Correlated 2n events')
plt.hist(data_cross, np.arange(20,135,binwidth), histtype='step', linewidth=1.5, label='Cross-talk events')
plt.legend()
# plt.grid()
plt.subplots_adjust(bottom=0.13)

plt.xlabel(r'$\theta{nn}$')
plt.ylabel(r'counts per bin')
plt.axes().set_axisbelow(True)
plt.show()