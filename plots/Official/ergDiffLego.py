import ROOT
import numpy as np
import mytools2 as mt2
from TH1Wrapper import TH2F, TH1F


#  ==================
n_erg_bins = 4
target = "DU"
theta_nbins = 12
# ==================
_f = ROOT.TFile("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/new_trees/DU.root")
tree_SP, n_pulses_SP = mt2.NCorrRun('SP', target, Forward=True, generate_dictionary=True).neutrons_doubles_tree
tree_DP, n_pulses_DP = mt2.NCorrRun('DP', target, Forward=True).neutrons_doubles_tree
tree_DP_multihit = _f.Get("tree")


erg_histSP = TH1F(0, 6, 100)
erg_histDP = TH1F(0, 6, 100)

erg_histDP.Project(tree_DP, "neutrons.coinc_hits[].erg", weight=1.0/n_pulses_DP)
erg_histSP.Project(tree_SP, "neutrons.coinc_hits[].erg", weight=1.0/n_pulses_SP)

erg_histSP += erg_histDP

bins, _ = mt2.median(erg_histSP, n_erg_bins)

print ("erg bins {}".format(bins))

histSPLego = TH2F(binarrays=bins)
histSPLego_w = TH2F(binarrays=bins)
histDPLego = TH2F(binarrays=bins)
histDPLegomultihit = TH2F(binarrays=bins)
histDPLego_w = TH2F(binarrays=bins)

histDPLego.Project(tree_DP, "neutrons.coinc_hits[].erg[0]:neutrons.coinc_hits[].erg[1]",weight=1.0/n_pulses_DP)
histSPLego.Project(tree_SP, "neutrons.coinc_hits[].erg[0]:neutrons.coinc_hits[].erg[1]", weight=1.0/n_pulses_SP)
histDPLegomultihit.Project(tree_DP_multihit, "erg[0]:erg[1]")

class DP_weights():
    def __init__(self):
        h_sp = histSPLego.__copy__()
        h_dp = histDPLegomultihit.__copy__()

        if target != "Cf252":
            h_sp -= 0.5*histDPLego

        h_sp.normalize(False)
        h_dp.normalize(False)

        self.w_h = h_sp/h_dp

        self.w_h.normalize()

    def get_weight(self, e1, e2):
        bin = self.w_h.FindBin(e1, e2)
        return self.w_h.GetBinContent(bin)


W = DP_weights()
# for evt in tree_DP:
#     if len(evt.neutrons.coinc_hits):
#         e1 = evt.neutrons.coinc_hits[0].erg[0]
#         e2 = evt.neutrons.coinc_hits[0].erg[1]
#         # print evt.neutrons.coinc_hits.nhits
#         # histDPLego.Fill(e1, e2)

if target != "Cf252":
    histcorrLego = histSPLego - 0.5*histDPLego
else:
    histcorrLego = histSPLego

histcorrLego.normalize()
histDPLegomultihit.normalize()
histDPLego.normalize()

# denom = histcorrLego.__copy__()

histcorrLego = histcorrLego/histDPLegomultihit
# histcorrLego /= denom

histcorrLego += histcorrLego.transpose()
histSPLego_w += histSPLego_w.transpose()

histcorrLego /= 2
histSPLego_w /= 2

# histcorrLego *= .8

histcorrLego.Draw("lego")
histcorrLego.GetZaxis().SetTitleOffset(1.2)
histcorrLego.SetLineWidth(2)
mt2.thesis_plot(histcorrLego, big_font=0.05)

histcorrLego.GetXaxis().SetTitle("E_{1} [MeV]")
histcorrLego.GetYaxis().SetTitle("E_{2} [MeV]")
histcorrLego.GetZaxis().SetTitle("nn_{corr}/nn_{uncorr}")
histcorrLego.SetStats(0)

histSP_theta = TH1F(20,180,theta_nbins)
histDP_theta = TH1F(20,180,theta_nbins)
histDP_w_theta = TH1F(20,180,theta_nbins)

histSP_theta.Project(tree_SP, "180/3.1415*neutrons.coinc_hits[].coinc_theta", weight=1.0/n_pulses_SP)
histDP_theta.Project(tree_DP, "180/3.1415*neutrons.coinc_hits[].coinc_theta", weight=1.0/n_pulses_DP)
histSP_w_theta = histSP_theta.__copy__()

for evt in tree_DP:
    for coinc_hit in evt.neutrons.coinc_hits:
        theta = 180/3.1415*coinc_hit.coinc_theta
        e1 = evt.neutrons.coinc_hits[0].erg[0]
        e2 = evt.neutrons.coinc_hits[0].erg[1]
        histDP_w_theta.Fill(theta, W.get_weight(e1, e2))

histDP_w_theta.update_bin_containers_from_hist()

histDP_w_theta *= sum(histDP_theta.binvalues)/sum(histDP_w_theta.binvalues)

# histDP_theta.Draw("hist ")
# histDP_w_theta.SetLineColor(ROOT.kRed)
# histDP_w_theta.Draw("same hist")

if target != "Cf252":
    histSP_theta -= 0.5*histDP_theta
histSP_theta /= histDP_theta

if target != "Cf252":
    histSP_w_theta -= 0.5*histDP_w_theta
histSP_w_theta /= histDP_w_theta

histSP_theta.normalize()
histSP_w_theta.normalize()
scale = 1.2/min(histSP_theta.binvalues[np.where(histSP_theta.binvalues>0)])
histSP_theta *= scale
histSP_w_theta *= scale

histSP_theta.SetLineColor(ROOT.kRed)

histSP_theta = histSP_theta.MySmooth(1)
histSP_w_theta = histSP_w_theta.MySmooth(1)

histSP_theta.binerrors *= 0.8
histSP_w_theta.binerrors *= 0.8

histSP_theta.Draw("E hist")
histSP_theta.SetMarkerStyle(33)
histSP_theta.SetLineStyle(2)
histSP_w_theta.SetMarkerStyle(24)
histSP_theta.GetXaxis().SetTitle("#theta_{nn} [degrees]")
histSP_theta.GetYaxis().SetTitle("nn_{corr}/nn_{uncorr}")
histSP_theta.SetMinimum(0)
histSP_theta.SetStats(0)
histSP_w_theta.Draw("same E hist", make_new_canvas=False)
mt2.thesis_plot(histSP_theta, big_font=0.06)

leg = ROOT.TLegend()
leg.AddEntry(histSP_theta, "unweighted", "lp")
leg.AddEntry(histSP_w_theta, "weighted", "lp")
leg.Draw()

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
