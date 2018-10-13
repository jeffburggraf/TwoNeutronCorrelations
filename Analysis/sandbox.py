import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt

binwidth = 12
# f = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/DP_FREYA_neutrons_coinc.root ")
# tree = f.Get("tree")

treeSP, n_pulsesSP = mt2.NCorrRun('SP',"DU").neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP',"DU").neutrons_doubles_tree


# def fuckyou():
#     global f1
#     f1 = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_DU/DP_weights.root")
#     weight_tree = f1.Get("tree")
#     print(weight_tree)
# # f1.Close()
#     return weight_tree
#
# weight_tree = fuckyou()
# print(weight_tree)
# # treeDP.AddFriend(weight_tree)


histSP = TH1F(0,10,binwidths=0.5)
drw = "0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])"
histSP.Project(treeSP, drw,  cut="neutrons.coinc_hits[0].ForwardDet[]==1")
# histDP.Project(treeDP, drw, cut="(1)", weight=1.0/n_pulsesDP)


# histSP -= 0.5*histDP
# histSP /= (histDP)

histSP.Draw()
#
#
#


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




