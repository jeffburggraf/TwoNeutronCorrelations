import ROOT
import numpy as np
import mytools2 as mt2
from TH1Wrapper import TH1F
import mytools as mt

target = "DU"
forward = True

treeSP, pulses_SP = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

treeDP, pulses_DP = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

n_bins = 20
histSP = TH1F(0,6,n_bins)
histDP = TH1F(0,6,n_bins)

drw = '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])'
# drw = 'abs(neutrons.coinc_hits[0].erg[0] -  neutrons.coinc_hits[0].erg[1])'
# drw = 'neutrons.coinc_hits[0].erg[]'
cut = '180/3.141*neutrons.coinc_hits.coinc_theta'
cut = mt.cut_rangeAND([0,180], cut)

if not forward:
    cut += ' && neutrons.coinc_hits.ForwardDet == 0'

print (histSP.Project(treeSP, drw, weight=1.0/pulses_SP, cut=cut))
histDP.Project(treeDP, drw, weight=1.0/pulses_DP, cut=cut)


# histSP -= 0.5*histDP

histDP.SetStats(0)
histSP.SetStats(0)
# ROOT.gStyle.SetOptStat('mn')

histSP.Draw()
histDP.Draw('sames')
histSP.GetXaxis().SetTitle('mean n-n energy [MeV]')
histSP.GetYaxis().SetTitle('Probability')

histDP.SetLineColor(ROOT.kRed)
histDP.SetMarkerStyle(20)
histSP.SetMarkerStyle(22)

histDP*=0.5

histSP *= 1.0/sum(histSP.binvalues)
histDP *= 1.0/sum(histDP.binvalues)

x_SP = histSP.bincenters[0][np.argmax(histSP.binvalues)]
x_DP = histDP.bincenters[0][np.argmax(histDP.binvalues)]


histSP.SetMaximum(1.15*max([max(histSP),max(histDP)]))

print(1/np.sqrt(2)*np.sqrt(np.sum([(np.sqrt(p) - np.sqrt(q))**2 for q,p in zip(histDP.binvalues, histSP.binvalues)])))


leg = ROOT.TLegend()

leg.AddEntry(histSP, 'Same pulse')
leg.AddEntry(histDP, 'Different pulse')
leg.Draw()

mt2.thesis_plot([histSP], big_font=True)

print (np.outer(histSP.binvalues, histDP.binvalues))
np.outer(histSP.binvalues, histDP.binvalues)


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
                        print(sys.exc_info())
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()



