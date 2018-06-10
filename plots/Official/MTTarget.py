import ROOT
import numpy as np
import matplotlib as mpl
import mytools as mt

mpl.use('TkAgg')

import mytools2 as mt2

from TH1Wrapper import TH1F


neutron_target = 'DU'
max_events = 3E6
forward = True

treeMT, _ = mt2.NCorrRun("SP", "MT", generate_dictionary=False, Forward=True).all_singles_tree
treeDU, _ = mt2.NCorrRun("SP", neutron_target, generate_dictionary=False, Forward=True).all_singles_tree
treeAl, pulses_Al = mt2.NCorrRun("SP", "Al", generate_dictionary=False, Forward=True).all_singles_tree

b = 1
histMT = TH1F(-10, 150, binwidths=b)
histDU = TH1F(-10, 150, binwidths=b)
histAl = TH1F(-10, 150, binwidths=b)

cut_MT = mt2.get_good_run_cut('MT')
cut_DU = mt2.get_good_run_cut(neutron_target)

print(cut_MT)
print(cut_DU)


DU_max = [1E6]

pulses_DU = mt2.get_pulses(cut_DU, neutron_target, max_events=DU_max)
pulses_MT = mt2.get_pulses(cut_MT, 'MT')

print(pulses_DU, DU_max)
print(pulses_MT)


c2 = "1"
if forward==False:
    c2 = "all.hits.ForwardDet == 0"

histMT.Project(treeMT, "all.hits.tof",weight=1.0/pulses_MT,cut=cut_MT + "&& " + c2)
histDU.Project(treeDU, "all.hits.tof",weight=1.0/pulses_DU,cut=cut_DU + "&& " + c2, max_events=DU_max[0])
histAl.Project(treeAl, "all.hits.tof",weight=1.0/pulses_Al,cut = c2)

histDU.SetLineWidth(4)
histAl.SetLineWidth(4)
histMT.SetFillColor(ROOT.kBlack)
histMT.SetFillStyle(3444)
histAl.SetLineStyle(2)


histMT.Draw('hist')
histAl.Draw('same  hist')
histMT.SetStats(0)
histMT.GetXaxis().SetTitle('ToF')
histMT.GetYaxis().SetTitle('counts per pulse')
# histDU.Draw('hist')
histMT.SetMaximum(max(histAl.binvalues)*1.15)

if neutron_target == 'D2O':
    neutron_target = 'D_{2}0'

leg1=ROOT.TLegend()
leg1.AddEntry(histMT, 'Empty', 'f')
leg1.AddEntry(histAl, 'Al', 'l')
leg1.Draw()

mt2.thesis_plot([histMT], big_font=1)


histDU.Draw('hist')
histAl.Draw('hist same')

histDU.GetXaxis().SetTitle('ToF')
histDU.GetYaxis().SetTitle('counts per pulse')

histDU.SetStats(0)
leg2=ROOT.TLegend()
leg2.AddEntry(histDU, neutron_target, 'f')
leg2.AddEntry(histAl, 'Al', 'f')
leg2.Draw()

mt2.thesis_plot([histDU], big_font=1)


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
