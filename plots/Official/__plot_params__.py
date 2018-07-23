

############## settings ##################

import matplotlib.pyplot as plt
import matplotlib as mpl
font = {'family':'DejaVu Sans',
        'size': 20}
mpl.rc('font', **font)
mpl.rc("savefig", dpi=400)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command

mpl.rc('text', usetex=True)


################ Error bar plot #####################
"""
plt.errorbar(x,y,err, label=label, marker=markers[index], markersize=6, fmt='o', capsize=3, elinewidth=.6,markeredgewidth=.6, c=colors[index])
"""


##################### Legend ##########################
"""
plt.legend((line,), (label,),title='Title', fontsize=13, bbox_to_anchor=(0.45, 1.02, 1.0-0.45, .5), loc=3, mode="expand", borderaxespad=0., frameon=False)
"""