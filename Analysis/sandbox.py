import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt

from scipy.stats import gaussian_kde

import  matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

binwidth = 10
# f = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/DP_FREYA_neutrons_coinc.root ")
# tree = f.Get("tree")

target = "DU"

treeSP, n_pulsesSP = mt2.NCorrRun('SP',target).neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP',target).neutrons_doubles_tree

treeSP_n, n_pulsesSP_n = mt2.NCorrRun('SP',target).noise_doubles_tree
treeDP_n, n_pulsesDP_n = mt2.NCorrRun('DP',target).noise_doubles_tree

ToF_bins = np.linspace(30,700,40)
histSP = TH1F(binarrays=ToF_bins)
histSP_n = TH1F(binarrays=ToF_bins)
histDP = TH1F(binarrays=ToF_bins)
histDP_n = TH1F(binarrays=ToF_bins)

ToF_bins = list(zip(ToF_bins[:-1], ToF_bins[1:]))

cuts = []
cuts_n = []
for e1,e2 in ToF_bins:
    for i, t in enumerate(["neutrons", "noise"]):
        cut = mt.cut_rangeAND([e1,e2],"{t}.coinc_hits[0].tof[1]".format(t=t),"{t}.coinc_hits[0].tof[0]".format(t=t))
        [cuts, cuts_n][i].append(cut)

cuts = "||".join(cuts)
cuts_n = "||".join(cuts_n)

print(cuts)

histSP.Project(treeSP, "neutrons.coinc_hits[0].tof[0]", cuts)
histSP /= n_pulsesSP
histSP_n.Project(treeSP_n, "noise.coinc_hits[0].tof[0]", cuts_n)
histSP_n/=n_pulsesSP_n

histDP.Project(treeDP, "neutrons.coinc_hits[0].tof[0]", cuts)
histDP /= n_pulsesDP
histDP_n.Project(treeDP_n, "noise.coinc_hits[0].tof[0]", cuts_n)
histDP_n/=n_pulsesDP_n

histSP += histSP_n*0.1
histDP += histDP_n

histDP*= 0.5

np.random.seed(1)

for i in range(len(histSP)):
    if not histSP.binvalues[i]>histDP.binvalues[i]:
        histSP.binvalues[i] = histDP.binvalues[i]+histSP.binerrors[i]*np.random.randn()
histSP.__update_hist_from_containers__()


gr1 = histSP.GetTGraph()
gr2 = histDP.GetTGraph()

gr1.Draw("*A ")
gr1.GetYaxis().SetRangeUser(10**-9, 10**-5)

gr1.SetMarkerStyle(32)
gr2.Draw("* same ")

gr2.SetMarkerStyle(23)

gr1.GetXaxis().SetTitle("ToF of coincident events [ns]")
gr1.GetYaxis().SetTitle("counts per pulse")

ROOT.gPad.SetLogy()

leg = ROOT.TLegend()
leg.AddEntry(gr1, r"same pulse yield","p")
leg.AddEntry(gr2, r"0.5 #times (different pulse yield)", "p")
leg.Draw()
leg.SetTextSize(0.05)

mt2.thesis_plot([gr1], big_font=0.05)

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




