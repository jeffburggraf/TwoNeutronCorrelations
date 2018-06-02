import ROOT
import numpy as np
import matplotlib as mpl
import mytools as mt

mpl.use('TkAgg')
import mytools2 as mt2

from TH1Wrapper import TH1F
from matplotlib import pyplot as plt
from matplotlib import pyplot as plt



font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size'   : 28}
## for Palatino and other serif fonts use:

# mpl.rc('font', **font)
# mpl.rc('text', usetex=True)
# mpl.rc("savefig", dpi=300)

target = "Cf252"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

histSP = TH1F(20,180,binwidths=7)
histDP = TH1F(20,180,binwidths=7)
tb = ROOT
cut = '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])>1.5 && neutrons.coinc_hits[0].ForwardTopBot==0'
print histSP.Project(treeSP_doubles, '180/3.1415*neutrons.coinc_hits.coinc_theta', cut)
print histDP.Project(treeDP_doubles, '180/3.1415*neutrons.coinc_hits.coinc_theta', cut, max_events=None)

hist_norm = (histSP)/histDP
hist_norm.SetStats(0)

ROOT.gStyle.SetOptStat('e');

ROOT.TGaxis.SetMaxDigits(3)

hist_norm.Draw('histE')
hist_norm.GetXaxis().SetTitle('#theta_{nn}')
hist_norm.GetYaxis().SetTitle('correlation [arb. units]')

hist_norm.SetLineWidth(2)

hist_norm.SetMinimum(0)

histSP.SetStats(0)
histSP.Draw('histE')
histSP.SetLineWidth(2)
histSP.GetXaxis().SetTitle('#theta_{nn}')
histSP.GetYaxis().SetTitle('counts')


mt2.thesis_plot([histSP,hist_norm], True, canvii=TH1F.tcanvii_refs)



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



