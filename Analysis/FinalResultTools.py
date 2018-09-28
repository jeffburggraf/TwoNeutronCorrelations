import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt


def get_weighted_cut(added_cut):
    treeSP, n_pulsesSP = mt2.NCorrRun('SP', "DU").neutrons_doubles_tree
    treeDP, n_pulsesDP = mt2.NCorrRun('DP', "DU").neutrons_doubles_tree

    erg_bins = np.linspace(0, 6, 5)

    global histSP, histDP
    histSP = TH2F(binarrays=erg_bins)
    histDP = TH2F(binarrays=erg_bins)

    histSP.Project(treeSP, "neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]", weight=1.0 / n_pulsesSP)
    histDP.Project(treeDP, "neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]", weight=1.0 / n_pulsesDP)

    histSP -= histDP * 0.5

    histSP += histSP.transpose()
    histSP /= 2.
    histSP.SetTitle("Corr")

    histDP += histDP.transpose()
    histDP /= 2.

    weight_hist = histSP.__copy__()
    weight_hist /= histDP

    cuts = []
    N = len(erg_bins) - 1

    for i in range(N):
        for j in range(i, N):
            erg1 = "neutrons.coinc_hits[0].erg[0]"
            erg2 = "neutrons.coinc_hits[0].erg[1]"

            _min1 = weight_hist.__binLeftEdges__[0][i]
            _max1 = weight_hist.__binRightEdges__[0][i]

            _min2 = weight_hist.__binLeftEdges__[1][j]
            _max2 = weight_hist.__binRightEdges__[1][j]

            cut = mt.cut_AND(mt.cut_rangeAND([_min1, _max1], erg1),
                             mt.cut_rangeAND([_min2, _max2], erg2))
            cut += " || " + mt.cut_AND(mt.cut_rangeAND([_min1, _max1], erg2),
                                       mt.cut_rangeAND([_min2, _max2], erg1))
            weight = round(weight_hist.binvalues[i][j], 2)
            cut = "{0}*({1})".format(weight, cut)
            cuts.append(cut)

    weightedDP_cut = "({1})*({0})".format("+".join(cuts), added_cut)

    if __name__ == "main":
        histSP.Draw("lego")
        histDP.Draw("lego")
        weight_hist.Draw('lego')

    return weightedDP_cut



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





