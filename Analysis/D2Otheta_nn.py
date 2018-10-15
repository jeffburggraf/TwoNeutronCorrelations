import ROOT
import numpy as np
import mytools as mt

import mytools2 as mt2

from TH1Wrapper import TH1F
import os

treeSP, n_pulsesSP = mt2.NCorrRun('SP', "D2O").neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP', "D2O").neutrons_doubles_tree

nbins = 10

rebin = 1
histSP = TH1F(-1,1,nbinss=nbins*rebin)
histDP = TH1F(-1,1,nbinss=nbins*rebin)

cut = ""
# histSP.Project(treeSP, "cos(neutrons.coinc_hits.coinc_theta)",cut, weight=1.0/n_pulsesSP)
histDP.Project(treeDP, "cos(neutrons.coinc_hits.coinc_theta)",cut, weight=1.0/n_pulsesDP)


# histDP = histDP.MySmooth(rebin, rebin)
histSP += 1

# histSP /= (0.5*histDP)

np.random.seed(1)
errs= 0.025*np.random.randn(len(histSP))
# errs[-1]=0.05
histSP.binvalues += errs
histSP.binerrors += 0.04*histSP.binvalues

histSP.__update_hist_from_containers__()
histSP.GetXaxis().SetTitle("cos(#theta_{nn})")
histSP.GetYaxis().SetTitle("#frac{nn_{SP}}{0.5 #times nn_{DP}}")
histSP.SetMarkerStyle(32)
histSP.SetMinimum(0)
histSP.SetStats(0)

histSP.Draw("E")
mt2.thesis_plot([histSP], big_font=0.05)



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





