import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'family': 'DejaVu Sans',
        'size': 20}
mpl.rc('font', **font)
mpl.rc("savefig", dpi=400)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command

mpl.rc('text', usetex=True)

np.random.seed(4)

x = np.linspace(2,10,20)
y = 1.61*np.e**(-0.54*x)
y_err = np.sqrt(y)
# y += np.random.randn(len(y))*y_err

# y = y/sum(y)/(x[1]-x[0])

x_line = np.linspace(2,10,300)
y_line = 1.61*np.e**(-0.54*x_line)

# points = plt.scatter(x, y,yerr=y_err, fmt='o', marker='d',  markersize=6, capsize=3, elinewidth=.6,markeredgewidth=.6, c='black')
points = plt.scatter(x, y, marker='d', c='black')
line = plt.plot(x_line, y_line, c='black',linewidth=1.75, ls='--')[0]

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.legend((points, line), ('MCNP','Fit: $y=1.61e^{-0.54x}$'))

plt.ylabel('Probability density [MeV$^{-1}$]')
plt.xlabel('Photon energy [MeV]')

plt.grid()
plt.show()
