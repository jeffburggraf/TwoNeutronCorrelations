import ROOT
import mytools as mt
import mytools2 as mt2
import numpy as np
from TH1Wrapper import TH1F

treeSP, pulses_SP = mt2.NCorrRun("SP","DU").neutrons_doubles_tree
treeDP, pulses_DP = mt2.NCorrRun("DP","DU").neutrons_doubles_tree

_f = ROOT.TFile("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/new_trees/DU.root")
treeDP_Dale = _f.Get("tree")

tb = ROOT.TBrowser()

nbins = 30
erg_hist_avg_SP = TH1F(0,10,nbins)
erg_hist_avg_DP_dale = TH1F(0,10,nbins)
erg_hist_avg_DP = TH1F(0,10,nbins)
erg_hist_diff_SP = TH1F(0,10,nbins)
erg_hist_diff_DP_dale = TH1F(-10,10,2*nbins)
erg_hist_diff_DP = TH1F(-10,10,2*nbins)

erg_hist_avg_SP.Project(treeSP, "0.5*(neutrons.coinc_hits[].erg[0] + neutrons.coinc_hits[].erg[1])")
erg_hist_avg_DP.Project(treeDP, "0.5*(neutrons.coinc_hits[].erg[0] + neutrons.coinc_hits[].erg[1])")
erg_hist_avg_DP_dale.Project(treeDP_Dale, "0.5*(erg[0] + erg[1])")

erg_hist_diff_SP.Project(treeSP, "abs(neutrons.coinc_hits[].erg[0] - neutrons.coinc_hits[].erg[1])")
erg_hist_diff_DP.Project(treeDP, "abs(neutrons.coinc_hits[].erg[0] - neutrons.coinc_hits[].erg[1])")
erg_hist_diff_DP_dale.Project(treeDP_Dale, "abs(erg[0] - erg[1])")


erg_hist_avg_SP.normalize()
erg_hist_avg_DP.normalize()
erg_hist_avg_DP_dale.normalize()

erg_hist_diff_SP.normalize()
erg_hist_diff_DP.normalize()
erg_hist_diff_DP_dale.normalize()

c1 = ROOT.TCanvas()
c1.Divide(2)
c1_1 = c1.cd(1)

ROOT.TGaxis.SetMaxDigits(3)


erg_hist_avg_SP.Draw("hist E",  make_new_canvas=False)

erg_hist_avg_SP.GetXaxis().SetTitle("mean neutron energy [MeV]")
erg_hist_avg_SP.GetYaxis().SetTitle("probability")

erg_hist_avg_DP.Draw("same E")
erg_hist_avg_DP_dale.Draw("same E ")
mt2.thesis_plot([erg_hist_avg_SP], big_font=0.05, Ytitle__offset=1.4)

c1.cd(2)
erg_hist_diff_SP.Draw("hist E",  make_new_canvas=False)

erg_hist_diff_SP.GetXaxis().SetTitle("absolute neutron energy difference [MeV]")
erg_hist_diff_SP.GetYaxis().SetTitle("probability")

erg_hist_diff_DP.Draw("same E")

erg_hist_diff_DP_dale.Draw("same E")

for (hist1, hist2, hist3) in [(erg_hist_avg_SP, erg_hist_avg_DP, erg_hist_avg_DP_dale), (erg_hist_diff_SP, erg_hist_diff_DP, erg_hist_diff_DP_dale)]:
    hist1.SetLineColor(ROOT.kRed)
    hist1.binerrors *= 0.8
    hist1.__update_hist_from_containers__()
    hist1.SetStats(0)
    hist1.SetLineWidth(3)
    hist2.SetMarkerStyle(27)
    hist3.SetMarkerStyle(20)

mt2.thesis_plot([erg_hist_diff_SP], big_font=0.05)


c2 = ROOT.TCanvas()
leg = ROOT.TLegend()
leg.AddEntry(erg_hist_avg_SP, "correlated neutron pairs (nn_{corr})", "l")
leg.AddEntry(erg_hist_avg_DP, "different pulse neutron pairs", "p")
leg.AddEntry(erg_hist_avg_DP_dale, "#splitline{different pulse neutron pairs}{from pulses with 2 neutron events (nn_{uncorr})}", "p")
leg.SetTextSize(0.05)
leg.Draw()



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

