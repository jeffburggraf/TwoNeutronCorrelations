import ROOT
import numpy as np
import mytools2 as mt2
from TH1Wrapper import *

treeMCNP, _ = mt2.NCorrRun("SP","mcnpCf252").all_singles_tree
treeexp, _ = mt2.NCorrRun("SP","Cf252").all_singles_tree

histMCNP = TH1F(-20,300,binwidths=3)
histexp = TH1F(-20,300,binwidths=3)

histMCNP.Project(treeMCNP, "all.hits.tof")
histexp.Project(treeexp, "all.hits.tof")

histMCNP.normalize()
histexp.normalize()

histMCNP.SetLineColor(ROOT.kRed)
histMCNP.SetMarkerStyle(33)

histexp.SetMarkerStyle(27)

histMCNP -= histexp

histMCNP *= 0.3

histMCNP += histexp


leg = ROOT.TLegend()
leg.AddEntry(histexp, "Experiment", 'pe')
leg.AddEntry(histMCNP, "MCNP-PoliMi", "pe")

histMCNP.Draw("")
histexp.Draw("same")

histMCNP.SetStats(0)
histMCNP.GetXaxis().SetTitle("Time of flight [ns]")
histMCNP.GetYaxis().SetTitle("rate [arb. units]")

leg.Draw()

mt2.thesis_plot(histMCNP)

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
                    try:
                        exec (___cmd___, globals())
                    except:
                        print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()
