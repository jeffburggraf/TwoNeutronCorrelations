import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt

import sys
sys.path.append('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis')
from FinalResultTools import*

# ====================
binwidths = [17, 16, 16, 13, 13]
n_erg_bins = 5
smooth = True
# ===================

rebin_factor = 4

binwidths = np.array(binwidths)
nbinss = list(map(np.int,(180-24)/binwidths))

_min_bin = 24

# if (180-_min_bin)%(binwidth) !=0:
#     _min_bin += + (180-_min_bin)%binwidth
#     print ("warning, raising lower bin to {}".format(_min_bin))
#     assert (180-_min_bin)%binwidth == 0

treeSP, n_pulsesSP = mt2.NCorrRun('SP',"DU").neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP',"DU").neutrons_doubles_tree

erg_hist = TH1F(0.4,6,100)
erg_hist.Project(treeSP, "neutrons.hits.erg")
erg_bins, _ = mt2.median(erg_hist, n_erg_bins)
del erg_hist

c1 = ROOT.TCanvas()
c1.Divide(n_erg_bins,2)

histos = []

_max = 0

for index, (E1, E2) in enumerate(zip(erg_bins[0:-1], erg_bins[1:])):
    c1.cd(index + 1)
    nbins = nbinss[index]*rebin_factor

    cut= mt.cut_rangeAND([E1,E2], "0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])")
    cut_DP= get_weighted_cut(cut)

    histSP = TH1F(_min_bin, 180, nbinss=nbins)
    histDP = TH1F(_min_bin, 180, nbinss=nbins)
    histDP_weighted = TH1F(_min_bin, 180, nbinss=nbins)

    drw = "180/3.1415*neutrons.coinc_hits[0].coinc_theta"
    histSP.Project(treeSP, drw, cut=cut, weight=1.0/n_pulsesSP)
    histDP.Project(treeDP, drw, cut=cut, weight=1.0/n_pulsesDP)
    histDP_weighted.Project(treeDP, drw, cut=cut_DP)
    histDP_weighted /= n_pulsesDP

    if smooth:
        histSP = histSP.MySmooth(rebin_factor, int(rebin_factor))
        histDP = histDP.MySmooth(rebin_factor, int(rebin_factor))
        histDP_weighted = histDP_weighted.MySmooth(rebin_factor, int(rebin_factor))

    histSP.SetMinimum(0)
    histDP.SetMinimum(0)

    histSP_old = histSP.__copy__()
    histSP_old -= 0.5*histDP
    histSP_old /= (0.5*histDP)

    histSP -= 0.5*histDP
    histSP /= (0.5*histDP_weighted)

    histSP.SetTitle("{0:.2f}<E<{1:.2f}".format(E1,E2))
    histSP_old.SetTitle("{0:.2f}<E<{1:.2f}".format(E1,E2))
    histSP.Draw(make_new_canvas= False)


    histos.append(histSP)
    histos.append(histSP_old)

    for hist in [histSP, histSP_old]:
        if max(hist.binvalues)>_max:
            _max = max(hist.binvalues)

    c1.cd(len(erg_bins)+index)

    histSP_old.SetMinimum(0)

    histSP_old.Draw(make_new_canvas=False)


for hist in histos:
    hist.SetMaximum(1.2*_max)



if __name__ == "__main__":
    import ROOT as ROOT
    from multiprocessing import Process, Queue
    import time, sys, os


    def input_thread(q, stdin):
        while True:
            print ('ROOT: ')
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
                    try:
                        exec (___cmd___, globals())
                    except:
                        print (sys.exc_info())
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()

