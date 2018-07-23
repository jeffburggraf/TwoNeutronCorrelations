import ROOT
import numpy as np
import mytools as mt

import mytools2 as mt2

from TH1Wrapper import TH1F

target = "DU"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

binwidth = 8
histSP = TH1F(20,180,binwidths=binwidth/3)
histDP = TH1F(20,180,binwidths=binwidth/3)
tb = ROOT
# cut = '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])>0.5'
cut = ''
print histSP.Project(treeSP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', cut, weight=1.0/pulses_SP_doubles)
print histDP.Project(treeDP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', cut, weight=1.0/pulses_DP_doubles)



hist_unnorm = histSP.Rebin(3)

if target == 'DU':
    hist_unnorm -= histDP.Rebin(3)*0.5

histSP.MySmooth(3.5)
histDP.MySmooth(3.5)

histSP = histSP.Rebin(3)
histDP = histDP.Rebin(3)

hist_norm = histSP.__copy__()
if target=='DU':
    hist_norm -= 0.5*histDP
hist_norm /= histDP

hist_norm.SetStats(0)

hist_norm *= 0.8/max(hist_norm.binvalues)


import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

plt.figure(figsize=(5, 10))


font = {'family':'DejaVu Sans',
        'size': 20}
mpl.rc('font', **font)
mpl.rc("savefig", dpi=400)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command

mpl.rc('text', usetex=True)
ax1 = plt.subplot(2,1,1)

plt.errorbar(hist_unnorm.x, hist_unnorm.binvalues, yerr=hist_unnorm.binerrors, linewidth=1, drawstyle='steps-mid', elinewidth=1., mec='black', capsize=2, c='black')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# plt.xlabel(r'$\theta _{nn}$')
plt.ylabel(r'$n$-$n$$_{\text{corr.}}$')

plt.xticks(np.arange(30, 180+30, 30))
y_ticks = list(map(lambda x: float('{:.0E}'.format(x)),np.linspace(0,max(hist_unnorm),5)))
plt.yticks(y_ticks)
plt.ylim(0,max(hist_unnorm.binvalues*1.15))
# plt.xlim(hist_unnorm.__binLeftEdges__[0][0], hist_unnorm.__binRightEdges__[0][-1])


plt.subplots_adjust(left=0.20)

plt.minorticks_on()

plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = plt.subplot(2,1,2, sharex=ax1)
ax2.grid(linestyle='-')
ax1.grid(linestyle='-')

plt.errorbar(hist_norm.bincenters[0], hist_norm.binvalues, yerr=hist_norm.binerrors,linewidth=1, drawstyle='steps-mid', elinewidth=1., mec='black', capsize=2, c='black')

plt.xlabel(r'$\theta _{nn}$')
plt.ylabel(r'$(n\text{-}n_{\text{corr.}})/(n\text{-}n_{\text{uncorr.}})$')
plt.minorticks_on()
plt.yticks(list(map(lambda x:mt2.round_to_n(x,1),np.linspace(0, mt2.round_to_n(max(hist_norm),1), 5))))
plt.ylim(0, max(hist_norm.binvalues*1.15))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# plt.xlim(hist_norm.__binLeftEdges__[0][0], hist_norm.__binRightEdges__[0][-1])

plt.savefig('/Users/jeffreyburggraf/PycharmProjects/2nCorrPhysRev/SPDPNormalization.png', transparent=True)
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



