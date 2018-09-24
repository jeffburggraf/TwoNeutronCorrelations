import numpy as np
from TH1Wrapper import TH1F,TH2F
import mytools2 as mt2
import os
from itertools import product
import wpca
import ROOT
import mytools as mt


bin_width = 15
tree_SP, n_pulses_SP = mt2.NCorrRun('SP','D2O', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP = mt2.NCorrRun('DP','D2O', Forward=True).neutrons_doubles_tree

histSP = TH1F(24,180,binwidths=bin_width/3)
histDP = TH1F(24,180,binwidths=bin_width/3)


drw = '180/3.1415*neutrons.coinc_hits.coinc_theta'
histSP.Project(tree_SP, drw, weight=1.0/n_pulses_SP)
histDP.Project(tree_DP, drw, weight=1.0/n_pulses_DP)

histSP.MySmooth(1)
histDP.MySmooth(1)

histSP = histSP.Rebin(3)
histDP = histDP.Rebin(3)

histSP /= (histDP*0.5)


histSP.SetMinimum(0)

histSP.GetYaxis().SetTitle('(n-n_{corr.})/(n-n_{uncorr.})')
histSP.GetXaxis().SetTitle('#theta_{nn}')

histSP.Draw()

mt2.thesis_plot([histSP], big_font=True)




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
                        print(sys.exc_info())
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()


#



