import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt


treeSP, n_pulsesSP = mt2.NCorrRun('SP',"DU").neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP',"DU").neutrons_doubles_tree

erg_hist = TH1F(0,6,150)
erg_hist.Project(treeSP, "neutrons.hits.erg")
erg_bins, _= mt2.median(erg_hist, 5)
erg_bins = np.linspace(0,6,5)
print(erg_bins)

histSP = TH2F(binarrays=erg_bins)
histDP = TH2F(binarrays=erg_bins)

histSP.Project(treeSP, "neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]", weight=1.0/n_pulsesSP)
histDP.Project(treeDP, "neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]",weight=1.0/n_pulsesDP)

histSP -= histDP*0.5



histSP += histSP.transpose()
histSP /= 2.
histSP.SetTitle("Corr")


histDP += histDP.transpose()
histDP /= 2.

histSP.SetMinimum(0)
histDP.SetMinimum(0)

histSP.Draw("lego")
histDP.Draw("lego")

W_hist = histSP/histDP

W_hist += W_hist.transpose()
W_hist *=0.5
W_hist.Draw("lego")

W_hist.SetMinimum(0)

cuts = []
N = len(erg_bins) - 1

weightedDP_cut = ""

for i in range(N):
    for j in range(i, N):
        erg1 = "neutrons.coinc_hits[0].erg[0]"
        erg2 = "neutrons.coinc_hits[0].erg[1]"
        _min1 = W_hist.__binLeftEdges__[0][i]
        _max1 = W_hist.__binRightEdges__[0][i]

        _min2 = W_hist.__binLeftEdges__[1][j]
        _max2 = W_hist.__binRightEdges__[1][j]

        cut = mt.cut_AND(mt.cut_rangeAND([_min1, _max1], erg1),
                         mt.cut_rangeAND([_min2, _max2], erg2))
        cut += " || "+mt.cut_AND(mt.cut_rangeAND([_min1, _max1], erg2),
                         mt.cut_rangeAND([_min2, _max2], erg1))
        weight = W_hist.binvalues[i][j]
        cut = "{0}*({1})".format(weight, cut)
        cuts.append(cut)
        # print(_min1,_max1, _min2, _max2, W_hist.binvalues[i][j])


weightedDP_cut = "({0})".format("+".join(cuts))

bin_width = 12
histSP_th = TH1F(24,180,binwidths=bin_width)
histDP_th = TH1F(24,180,binwidths=bin_width)
histDP_th_weighted = TH1F(24, 180,binwidths=bin_width)

drw = "180/3.1415*neutrons.coinc_hits.coinc_theta"
histSP_th.Project(treeSP, drw, weight=1.0/n_pulsesSP)
histDP_th.Project(treeDP, drw, weight=1.0/n_pulsesDP)
histDP_th_weighted.Project(treeDP, drw, cut=weightedDP_cut)
histDP_th_weighted /= n_pulsesDP

old_metheod_hist = histSP_th.__copy__()

histSP_th -= 0.5*histDP_th

histSP_th /= (0.5*histDP_th_weighted)
histSP_th.Draw("hist E")
histSP_th.SetMinimum(0)








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
