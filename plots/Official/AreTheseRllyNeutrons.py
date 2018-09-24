import ROOT
import numpy as np
import mytools as mt

import mytools2 as mt2

from TH1Wrapper import TH1F


#############################
binwidth = 1
#############################

def get_n_ph_combined_hist(target, min_max_binwidth, max_photon_events = 10E6, SPDP = 'SP', drw = '{par}.coinc_hits[].tof', cut = ''):
    assert '{par}' in drw

    tree_n, pulses_n = mt2.NCorrRun('DP' if target!='Mt' else 'SP', target).neutrons_doubles_tree
    tree_ph, _ = mt2.NCorrRun(SPDP, target).photons_doubles_tree

    if max_photon_events is not None:
        max_photon_events = int(max_photon_events)
    else:
        max_photon_events = tree_ph.GetEntries()

    _min, _max, binwidth = min_max_binwidth
    hist_n = TH1F(_min, _max, binwidths=binwidth)
    hist_ph = TH1F(_min, _max, binwidths=binwidth)

    tree_ph.GetEntry(max_photon_events)
    pulses_ph = tree_ph.PulseNumber

    hist_n.Project(tree_n, drw.format(par='neutrons'), cut=cut, weight=1.0/pulses_n)
    hist_ph.Project(tree_ph, drw.format(par='photons'), cut=cut, weight=1.0/pulses_ph, max_events=max_photon_events)

    for i in range(len(hist_ph)):
        if hist_ph.binvalues[-i] !=0:
            if hist_ph.binvalues[-(i+1)] !=0:
                b = hist_ph.binvalues[-(i+1)]
                I_ph = len(hist_ph)-(i+1)
                break
    else:
        b = None

    for i in range(len(hist_n)):
        if hist_n.binvalues[i] !=0:
            if hist_n.binvalues[i+1] != 0:
                a = hist_n.binvalues[i+1]
                I_n=i+1
                break
    else:
        a = None

    if None not in [a, b]:
        hist_ph *= a/b
        hist_ph.binerrors *= hist_n.binerrors[I_n]/hist_ph.binerrors[I_ph]
        hist_ph.__update_hist_from_containers__()

    hist_n += hist_ph

    return hist_n



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







