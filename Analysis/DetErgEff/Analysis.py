import ROOT
import numpy as np
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F
import os


treeSP, n_pulsesSP = mt2.NCorrRun('SP', "Cf252").neutrons_doubles_tree

nbins = 15


_max = 7
histSP = TH1F(0,_max,nbinss=nbins)
hist_theory = TH1F(0,_max,nbinss=nbins)

root_f=  ROOT.TFile("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/DetErgEff/ptrac0.root")
tree = root_f.Get("tree")
hist_theory.Project(tree, "erg","evt_type == 5")
cut = "neutrons.hits.ForwardDet == 0"
histSP.Project(treeSP, "neutrons.hits.erg",cut)



histSP /= hist_theory

func = ROOT.TF1("func","((x**2-[1])*2.718**(-[0]*x**2)+[1])")

func.SetParameter(0,0.38)
func.SetParameter(1,1.16)
#
histSP.Draw()
histSP.Fit("func")
histSP.SetMarkerStyle(33)


histSP.SetStats(0)
histSP.GetXaxis().SetTitle("Neutron Energy [MeV]")
histSP.GetYaxis().SetTitle("rel. efficiency [arb. units]")
mt2.thesis_plot([histSP], big_font=0.05)

hist_theory.Draw()
histSP.Draw()
leg = ROOT.TLegend()
leg.AddEntry(histSP, "measurement", "lp")
leg.AddEntry(func, "Fit: y = (x^{{2}} - {1:.1f})e^{{-{0:.2f}x^{{2}}}}+{1:.1f}".format(func.GetParameter(0), func.GetParameter(1)), "lp")
leg.Draw()
p = [leg]
histos = []


nbins = 7

hist_theory = TH1F(0,_max,nbins)
hist_theory.Project(tree, "erg","evt_type == 5")

legs= []
c1 = ROOT.TCanvas()
c1.Divide(2,2)
_angles = mt.angles[1:-1]
line_styles = [1,2,9,5]
for Index, angles in enumerate(np.array_split(_angles, 4)):
    leg = ROOT.TLegend()
    c1.cd(Index + 1)
    legs.append(leg)
    _max_value = 0
    grs = []
    for index, det in enumerate(angles):
        pad = c1.cd(Index + 1)
        hist = TH1F(0, _max, nbins)
        hist.Project(treeSP, "neutrons.hits.erg", "neutrons.hits.det == {}".format(det))
        hist /= hist_theory
        hist.normalize(False)
        hist *= 7
        gr = ROOT.TGraphErrors(
            len(hist.bincenters[0]), np.array(hist.bincenters[0],dtype=np.float64), np.array(hist.binvalues, dtype=np.float64),
        np.zeros_like(hist.binvalues, dtype=np.float64), hist.binerrors)

        gr.SetMarkerStyle(24 + index)
        histos.append(gr)
        gr.SetLineStyle(line_styles[index])
        leg.AddEntry(gr, str(det), "lp")
        leg.SetTextSize(0.1)
        grs.append(gr)
        gr.SetTitle("")
        if max(hist.binvalues)>_max_value:
            _max_value = max(hist.binvalues)
            gr.GetYaxis().SetRangeUser(0, 1.15 * _max_value)

    grs[0].Draw("Alp")
    mt2.thesis_plot([grs[0]], big_font=0.1)
    grs[0].GetXaxis().SetTitle("Neutron Energy [MeV]")
    grs[0].GetYaxis().SetTitle("rel. efficiency [arb. units]")
    grs[0].GetXaxis().SetTitleOffset(1)

    leg.Draw()
    for gr in grs[1:]:
        gr.Draw("lpsame")

        gr.SetMinimum(0)



l = ROOT.TLatex()
# l.DrawLatex(2,0.4,"y = (x^{{2}} - {1:.1f})e^{{-{0:.2f}x^{{2}}}}+{1:.1f}".format(func.GetParameter(0), func.GetParameter(1)));

s = []
for x,y in zip(histSP.bincenters[0], histSP.binvalues):
    s .append("{{ {0:.2f}, {1:.2f} }}".format(x,y))
print(" ,".join(s))

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






