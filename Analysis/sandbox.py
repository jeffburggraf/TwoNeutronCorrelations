import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt

from scipy.stats import gaussian_kde

import  matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

binwidth = 12
# f = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/DP_FREYA_neutrons_coinc.root ")
# tree = f.Get("tree")

treeSP, n_pulsesSP = mt2.NCorrRun('SP',"DU").noise_singles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP',"DU").noise_singles_tree

histSP = TH1F(24,700,binwidths=0.05)
histDP = TH1F(24,700,binwidths=0.05)

histSP.Project(treeSP,"noise.hits.tof")
histSP.Draw()






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




