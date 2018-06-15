
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

plt.figure(figsize=(10,20))


font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size': 30} # 18 for single figure, 30 for double
## for Palatino and other serif fonts use:

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc("savefig", dpi=400)

ax = plt.subplot(2,1,1)

plt.errorbar(histSP.bincenters[0], histSP.binvalues, yerr=histSP.binerrors, drawstyle='steps-mid', elinewidth=1.5, mec='black', capsize=4, c='black')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'$\theta _{nn}$')
plt.ylabel('counts/trigger')

plt.ylim(0,max(histSP.binvalues*1.15))
plt.minorticks_on()
plt.xticks(np.arange(0, 200, 30))
plt.legend(loc='upper left', fontsize=25)


plt.savefig('/Users/jeffreyburggraf/PycharmProjects/2nCorrPhysRev/test.png', transparent=True)
plt.show()