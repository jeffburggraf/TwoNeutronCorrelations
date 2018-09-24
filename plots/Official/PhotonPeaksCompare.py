import ROOT
import numpy as np
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F

target = "Al"

treeAl, pulses_Al = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).photons_singles_tree
treeDU, pulses_DU = mt2.NCorrRun("SP", 'DU', generate_dictionary=False, Forward=True).photons_singles_tree

histAL = TH1F(-20,40,binwidths=0.5)
histDU = TH1F(-20,40,binwidths=0.5)

histAL.GetXaxis().SetTitle('ToF [ns]')
histAL.GetYaxis().SetTitle('counts per pulse')
histAL.SetMarkerStyle(32)


AL_events = histAL.Project(treeAl, 'photons.hits.tof', weight=1.0/pulses_Al)
DU_events = histDU.Project(treeDU, 'photons.hits.tof', weight=1.0/pulses_DU)

scale = DU_events/pulses_DU/(AL_events/pulses_Al)
print('{1} scale is: {0:.2f}'.format(scale,target))

histAL *= scale


histAL.Draw('hist')
histAL.SetLineColor(ROOT.kRed)
histDU.Draw('hist same')

leg = ROOT.TLegend()
leg.AddEntry(histAL, target)
leg.AddEntry(histDU, 'DU')
leg.Draw()

mt2.thesis_plot([histAL])

treeAl, pulses_Al = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).photons_doubles_tree
treeDU, pulses_DU = mt2.NCorrRun("SP", 'DU', generate_dictionary=False, Forward=True).photons_doubles_tree

histAL_th = TH1F(0, 180, binwidths=1)
histDU_th = TH1F(0, 180, binwidths=1)

histAL_th.GetXaxis().SetTitle('#theta_{abs}')
histAL_th.GetYaxis().SetTitle('counts per pulse')

histAL_th.Project(treeAl, '180/3.1415*photons.coinc_hits.theta_abs', weight=1.0/pulses_Al)
histDU_th.Project(treeDU,  '180/3.1415*photons.coinc_hits.theta_abs', weight=1.0/pulses_DU)

histAL_th *= scale


histAL_th.Draw('hist')
histAL_th.SetLineColor(ROOT.kRed)
histDU_th.Draw('same hist')

leg = ROOT.TLegend()
leg.AddEntry(histAL_th, target)
leg.AddEntry(histDU_th, 'DU')
leg.Draw()

mt2.thesis_plot([histAL_th])



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






