import numpy as np
from TH1Wrapper import TH1F,TH2F
import mytools2 as mt2
import os
from itertools import product

import ROOT
import mytools as mt
import mytools2 as mt2

bin_width = 2

tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree

hist = TH1F(24,180,binwidths=bin_width)

hist.Project(tree_SP, "180/3.1415*neutrons.coinc_hits.coinc_theta",weight=1.0/n_pulses_SP/bin_width)


hist.Draw("hist E")
hist.SetLineWidth(2)
hist.GetXaxis().SetTitle("#theta_{nn}  [degrees]")
hist.GetYaxis().SetTitle("counts/pulse/degree")
mt2.thesis_plot(hist, big_font=0.055, Ytitle__offset=1)

hist.SetStats(0)

tb = ROOT.TBrowser()


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
