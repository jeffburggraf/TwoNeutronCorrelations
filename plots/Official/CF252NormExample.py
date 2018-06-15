import ROOT
import numpy as np
import mytools as mt

import mytools2 as mt2

from TH1Wrapper import TH1F


target = "Cf252"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

histSP = TH1F(20,180,binwidths=7)
histDP = TH1F(20,180,binwidths=7)
tb = ROOT
cut = '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])>1.5 && neutrons.coinc_hits[0].ForwardTopBot==0'
print histSP.Project(treeSP_doubles, '180/3.1415*neutrons.coinc_hits.coinc_theta', cut, weight=1.0/pulses_SP_doubles)
print histDP.Project(treeDP_doubles, '180/3.1415*neutrons.coinc_hits.coinc_theta', cut)

hist_norm = (histSP)/histDP
hist_norm.SetStats(0)

hist_norm *= 0.8/max(hist_norm.binvalues)


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

plt.xticks(np.arange(0, 200, 30))
plt.ylim(0,max(histSP.binvalues*1.15))
plt.xlim(histSP.__binLeftEdges__[0][0], histSP.__binRightEdges__[0][-1])


plt.subplots_adjust(bottom=0.17)

plt.minorticks_on()

ax = plt.subplot(2,1,2)

plt.errorbar(hist_norm.bincenters[0], hist_norm.binvalues, yerr=hist_norm.binerrors, drawstyle='steps-mid', elinewidth=1.5, mec='black', capsize=4, c='black')

plt.xlabel(r'$\theta _{nn}$')
plt.ylabel('correlation [arb. units]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.minorticks_on()
plt.xticks(np.arange(0, 200, 30))
plt.ylim(0, max(hist_norm.binvalues*1.15))
plt.xlim(hist_norm.__binLeftEdges__[0][0], hist_norm.__binRightEdges__[0][-1])

plt.savefig('/Users/jeffreyburggraf/PycharmProjects/2nCorrPhysRev/Cf252Norm.png', transparent=True)
plt.show()

if __name__ == "__main__":
    import ROOT as ROOT
    from multiprocessing import Process, Queue
    import time, sys, os
    def input_thread(q, stdin):
        while True:
            print 'ROOT: '
            cmd = stdin.readline()
            q.put(cmd)
    def root(char):
        assert isinstance(char, str), "Argument must be string!"
        ROOT.gROOT.ProcessLine(char)
    if __name__ == '__main__':
        ___queue___ = Queue()
        ___newstdin___ = os.fdopen(os.dup(sys.stdin.fileno()))
        ___input_p___ = Process(target=input_thread,
                                args=(___queue___, ___newstdin___))
        ___input_p___.daemon = True
        ___input_p___.start()
        ___g___ = ROOT.gSystem.ProcessEvents
        try:
            while 1:
                if not ___queue___.empty():
                    ___cmd___ = ___queue___.get()
                    if ___cmd___ == 'tb':
                        ___tb___ = ROOT.TBrowser()
                    else:
                        try:
                            exec (___cmd___, globals())
                        except:
                            print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()



