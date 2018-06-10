import ROOT
import numpy as np
import matplotlib as mpl
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F

mpl.use('TkAgg')

neutron_target = 'DU'

histDU_n = TH1F(-20,150,binwidths=1.5)
histDU_ph = TH1F(-20,150,binwidths=1.5)
histAl = TH1F(-20,150,binwidths=1.5)

treeDU_n, n_pulses = mt2.NCorrRun("SP", neutron_target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeDU_a, _ = mt2.NCorrRun("SP", neutron_target, generate_dictionary=False, Forward=True).photons_singles_tree
treeAl, pulses_Al = mt2.NCorrRun("SP", "Al", generate_dictionary=False, Forward=True).all_singles_tree

histDU_n.Project(treeDU_n, "neutrons.hits.tof", weight=1.0/n_pulses)
histDU_ph.Project(treeDU_a, "photons.hits.tof", max_events=1E6)
histAl.Project(treeAl, "all.hits.tof", max_events=1E6)

treeDU_a.GetEntry(int(1E6))
n_pulses_ph=treeDU_a.PulseNumber
histDU_ph *= 1.0/n_pulses_ph


s = histDU_ph.binvalues[histDU_ph.FindBin(25)]/ histDU_n.binvalues[histDU_n.FindBin(30)]

histDU_n *= s
histAl *= s

histAl *= 0.85*max(histDU_ph.binvalues)/max(histAl.binvalues)
histDU_n.SetMaximum(max(histDU_ph.binvalues)*1.15)

histDU_n += histDU_ph

histAl.SetFillStyle(3003)
histAl.SetFillColor(ROOT.kBlack)
histDU_n.SetLineWidth(3)
histDU_n.SetStats(0)

histDU_n.GetXaxis().SetTitle('ToF [ns]')
histDU_n.GetYaxis().SetTitle('counts/pulse')

histDU_n.Draw('hist')
histAl.Draw('hist same')

leg = ROOT.TLegend()
leg.AddEntry(histDU_n, 'DU','fl')
leg.AddEntry(histAl, 'Al', 'f')
leg.Draw()

ROOT.gPad.SetLogy()
mt2.thesis_plot([histDU_n], big_font=True)

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
