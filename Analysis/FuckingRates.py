import ROOT
import numpy as np
import mytools as mt

import mytools2 as mt2

from TH1Wrapper import TH1F
import os

from collections import OrderedDict


tree_neutrons_singles, n_pulses_neutrons_singles = mt2.NCorrRun('SP', "DU").neutrons_singles_tree
tree_neutrons_doubles, n_pulses_neutrons_doubles = mt2.NCorrRun('SP', "DU").neutrons_doubles_tree
tree_photons_singles, n_pulses_photons_singles = mt2.NCorrRun('SP', "DU").photons_singles_tree
tree_photons_doubles, n_pulses_photons_doubles = mt2.NCorrRun('SP', "DU").photons_doubles_tree
print(tree_neutrons_singles.Draw("","neutrons.ncoinc==1"), "dmfkf")
tb = ROOT.TBrowser()
singles_neutrons = OrderedDict({})
for det in mt.angles:
    topbots = [0] if det not in [30, 330] else [-1,1]
    for topbot in topbots:
        cut ="neutrons.hits.det == {0}{1}".format(det, "" if topbot == 0 else "&& neutrons.hits.ForwardDet == {}".format(topbot))
        cut_ph ="photons.hits.det == {0}{1}".format(det, "" if topbot == 0 else "&& photons.hits.ForwardDet == {}".format(topbot))
        print(det, cut)
        n_s_rate = float(tree_neutrons_singles.Draw("", cut))
        n_s_rate += tree_neutrons_doubles.Draw("",cut)
        ph_s_rate = tree_photons_singles.Draw("", cut_ph)
        ph_s_rate += n_pulses_photons_singles*tree_photons_doubles.Draw("",cut_ph)/n_pulses_photons_doubles
        if det in [30,330]:
            det_label = str(det) + (" top" if topbot == 1 else " bottom")
        else:
            det_label = str(det)
        singles_neutrons[det_label] = ("{:.2E}".format(n_s_rate/n_pulses_neutrons_singles),
                                       "{:.2E}".format(ph_s_rate/n_pulses_photons_singles))
    print("det ",det,"",  n_s_rate, " ", ph_s_rate)

s = []
for det, (rate_n,rate_ph) in singles_neutrons.iteritems():
    s.append("{0} & {1} & {2} \\\\".format(det,rate_n, rate_ph))

s = " \midrule\n".join(s)

print("\n\n\n ")
print (s)
print("\n")

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