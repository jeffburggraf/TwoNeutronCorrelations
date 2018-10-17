import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt

from scipy.stats import gaussian_kde

import  matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

binwidth = 12
# f = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/DP_FREYA_neutrons_coinc.root ")
# tree = f.Get("tree")

treeSP_noise, n_pulses_noiseSP = mt2.NCorrRun('SP',"DU").noise_doubles_tree
treeDP_noise, n_pulses_noiseDP = mt2.NCorrRun('DP',"DU").noise_doubles_tree
treeSP_neutons, n_pulses_neutronsSP = mt2.NCorrRun('SP',"DU").neutrons_doubles_tree
treeDP_neutons, n_pulses_neutronsDP = mt2.NCorrRun('DP',"DU").neutrons_doubles_tree

bw = 5
hist_neutronsSP = TH1F(30,700,binwidths=bw)
hist_neutronsDP = TH1F(30,700,binwidths=bw)
hist_noiseSP = TH1F(30,700,binwidths=bw)
hist_noiseDP = TH1F(30,700,binwidths=bw)


drw = "0.5*({0}.coinc_hits[0].tof[1] + {0}.coinc_hits[0].tof[0])"
hist_neutronsSP.Project(treeSP_neutons, drw.format("neutrons") , weight=1.0/n_pulses_neutronsSP)
hist_neutronsDP.Project(treeDP_neutons, drw.format("neutrons"), weight=1.0/n_pulses_neutronsDP)
hist_noiseSP.Project(treeSP_noise, drw.format("noise"), weight=1.0/n_pulses_noiseSP)
hist_noiseDP.Project(treeDP_noise, drw.format("noise"), weight=1.0/n_pulses_noiseDP)


hist_neutronsSP +=hist_noiseSP
hist_neutronsDP +=hist_noiseDP


hist_neutronsSP_raw = hist_neutronsSP.__copy__()
hist_neutronsSP_raw += hist_neutronsSP

# hist_neutronsSP -= 0.5*(hist_neutronsDP)
# hist_noiseSP -= 0.5*(hist_noiseDP)
#
# hist_neutronsSP += hist_noiseSP
hist_neutronsDP *=0.5

hist_neutronsSP_raw.Draw("hist")
hist_neutronsDP.Draw("hist same")
hist_neutronsDP.SetLineColor(ROOT.kRed)

hist_noiseSP.Draw()



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




