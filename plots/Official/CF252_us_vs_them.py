import ROOT
import numpy as np
import matplotlib as mpl
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size'   : 18}
## for Palatino and other serif fonts use:

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc("savefig", dpi=400)

neutron_target = 'Cf252'

treeSP, n_pulsesSP = mt2.NCorrRun("SP", neutron_target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun("DP", neutron_target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

histSP = TH1F(0, 180, binwidths=9)
histDP = TH1F(0, 180, binwidths=9)

histSP.Project(treeSP, '180/3.14*neutrons.coinc_hits[].coinc_theta',weight=1.0/n_pulsesSP)#, max_events=5E5)
histDP.Project(treeDP, '180/3.14*neutrons.coinc_hits[].coinc_theta', weight=1.0/n_pulsesDP)#, max_events=5E5)

histSP /= histDP

histSP *= 1.6 / (histSP.binvalues[-1])
histSP.SetMinimum(0.775)
histSP.SetMaximum(2.5)


plt.errorbar(histSP.bincenters[0], histSP.binvalues, yerr=histSP.binerrors, marker='s', linewidth = 0, elinewidth=1, mec='black', ms=6, mfc='blue', capsize = 3)

plt.ylim(0.775,2.5)

plt.ylabel(r'$\theta _{nn}$')
plt.xlabel('correlation [arb. units]')
plt.subplots_adjust(bottom=0.17)
# plt.grid()
plt.minorticks_on()
plt.xticks(np.arange(0, 200, 20))

# plt.locator_params(axis='x', nbins=6)
plt.savefig('test.png', transparent=True)

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
