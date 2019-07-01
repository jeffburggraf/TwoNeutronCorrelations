import ROOT
import numpy as np
import matplotlib as mpl
import mytools as mt

mpl.use('TkAgg')
import mytools2 as mt2

from TH1Wrapper import TH1F
from matplotlib import pyplot as plt
from matplotlib import pyplot as plt



font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size'   : 28}
## for Palatino and other serif fonts use:

# mpl.rc('font', **font)
# mpl.rc('text', usetex=True)
# mpl.rc("savefig", dpi=300)

target = "DU"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).all_singles_tree

histSP = TH1F(-10,150,binwidths=1.5)

print histSP.Project(treeSP_doubles, 'all.coinc_hits[].tof','all.coinc_hits[].ForwardDet==0', max_events=2E6, weight=1.0/2E6)

hist_photons = histSP.__copy__()
hist_neutrons = histSP.__copy__()

hist_photons.binvalues = np.where(np.array(hist_photons.bincenters[0])<25, hist_photons.binvalues, 0)
hist_photons.binvalues = np.where(np.array(hist_photons.bincenters[0])>-6, hist_photons.binvalues, 0)

hist_neutrons.binvalues = np.where(np.array(hist_neutrons.bincenters[0])>41, hist_neutrons.binvalues, 0)
hist_neutrons.binvalues = np.where(np.array(hist_neutrons.bincenters[0])<130, hist_neutrons.binvalues, 0)

hist_photons.__update_hist_from_containers__()
hist_neutrons.__update_hist_from_containers__()

ROOT.gStyle.SetOptStat('e');

ROOT.TGaxis.SetMaxDigits(3)

histSP.SetStats(0)
histSP.SetLineWidth(2)
histSP.GetXaxis().SetTitle('ToF')
histSP.GetYaxis().SetTitle('counts per pulse')


hist_photons.SetFillStyle(3244)
hist_neutrons.SetLineColorAlpha(ROOT.kWhite, 0)
hist_photons.SetLineColorAlpha(ROOT.kWhite, 0)

hist_neutrons.SetFillStyle(3350)
hist_photons.SetFillColor(ROOT.kBlack)
hist_neutrons.SetFillColor(ROOT.kBlack)

histSP.Draw('histE')
hist_photons.Draw('hist same')
hist_neutrons.Draw('hist same')
ROOT.gPad.SetLogy()
mt2.thesis_plot([histSP], True, canvii=TH1F.tcanvii_refs)

leg = ROOT.TLegend()
leg.AddEntry(hist_photons, 'Gammas', 'f')
leg.AddEntry(hist_neutrons, 'Neutrons', 'f')
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



