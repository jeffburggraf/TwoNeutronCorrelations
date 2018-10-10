import ROOT
import os
from TH1Wrapper import TH2F, TH1F
import re
import mytools2 as mt2
import time
import numpy as np

treeDP, _ = mt2.NCorrRun("DP","DU").neutrons_doubles_tree

f1 = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_DU/DP_weights.root")
weight_tree = f1.Get("tree")

# treeDP.AddFriend(weight_tree)

hist1 = TH1F(0,10,50)
hist2 = TH1F(0,10,50)

hist1.Project(treeDP, "neutrons.coinc_hits[].erg[]","DP_weight")
hist2.Project(treeDP, "neutrons.coinc_hits[].erg[]","")

hist1.Draw("hist E")
hist2.Draw("same hist E" )
hist2.SetLineColor(ROOT.kRed)

print(sum(hist2.binvalues))
print(sum(hist1.binvalues))

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


