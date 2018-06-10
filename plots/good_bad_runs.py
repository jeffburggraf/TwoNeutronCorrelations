import ROOT
import numpy as np
import matplotlib as mpl
import mytools as mt
import re

mpl.use('TkAgg')
import mytools2 as mt2

from TH1Wrapper import TH1F

target = 'D2O'

tree, pulses = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).all_singles_tree

runs_pulses = {}

for runnum in mt2.runs[target]:
    histMT = TH1F(-20, 150, binwidths=1, title=runnum)
    cut = "RunNumber == {}".format(runnum)
    histMT.Project(tree, "all.hits.tof",cut = cut, weight=1.0/pulses)

    histMT.SetFillColor(ROOT.kBlack)
    histMT.SetFillStyle(3144)

    histMT.Draw('hist')
    l = ROOT.TLegend()
    l.AddEntry(histMT, str(runnum), 'f')
    l.Draw()
    ROOT.gPad.SetLogy()

    raw_Tfile = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Aug_2017/r{}.root".format(runnum))
    _tree_ = raw_Tfile.Get("FF")

    runs_pulses[runnum] = int(_tree_.GetEntries())
    raw_Tfile.Close()

print("'{0}':{1}".format(target, runs_pulses))


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
                    if ___cmd___ == 'tb':
                        ___tb___ = ROOT.TBrowser()
                    else:
                        try:
                            exec (___cmd___, globals())
                        except:
                            print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()
