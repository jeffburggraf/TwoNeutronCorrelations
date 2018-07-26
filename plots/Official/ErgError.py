import matplotlib as mpl
mpl.use('TkAgg')

from matplotlib import pyplot as plt
import numpy as np

font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size'   : 28}
## for Palatino and other serif fonts use:

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc("savefig", dpi=300)

tofs = np.arange(35,145, 0.01)

erg = 8127/tofs**2
err = 40000/tofs**3

plt.figure(1, figsize=(10,8))
l2 = plt.plot(tofs, erg, linewidth=4,c='black')
plt.grid()
plt.xlabel('ToF [ns]')
plt.ylabel('Energy [MeV]')
# plt.subplots_adjust(bottom = 0.15)
plt.savefig('/Users/jeffreyburggraf/Pictures/ToF2Erg.png',bbox_inches='tight')



# plt.legend()
# plt.plot(tofs, erg+err, linestyle ='dashed', color = 'black')
# plt.plot(tofs, erg-err,linestyle ='dashed',color = 'black')
# coll = plt.fill_between(tofs,erg+err, erg-err, color = '#539caf', alpha = 0.2, linestyle ='-', label=r'erg$\pm \sigma$', lw =1, edgecolor = 'black')
# coll.set_edgecolor('black')
plt.figure(2, figsize=(10,8))
# plt.subplots_adjust(bottom = 0.15)

l1 = plt.plot(erg, err, linewidth=4, c='black')


plt.ylabel('$\Delta$Energy [MeV]')
plt.xlabel('Energy [MeV]')
plt.grid()
plt.locator_params(axis='x', nbins=10)

# plt.tight_layout()
plt.savefig('/Users/jeffreyburggraf/Pictures/DeltaErg.png', bbox_inches='tight')



plt.show()

