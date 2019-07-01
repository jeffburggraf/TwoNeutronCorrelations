import ROOT
import mytools2 as mt2
import numpy as np
from TH1Wrapper import TH1F


bin_width = 1

tree_n, n_pulses_n = mt2.NCorrRun("SP", "CF252").neutrons_singles_tree
tree_ph, n_pulses_ph = mt2.NCorrRun("SP", "CF252").photons_doubles_tree

hist_n = TH1F(-20, 120, binwidths=bin_width)
hist_ph = TH1F(-20, 29, binwidths=bin_width)
hist_fuck = TH1F(-20, 29, binwidths=bin_width)
hist_fuck.FillRandom("gaus", 10000)
hist_fuck.update_bin_containers_from_hist()
hist_fuck.binvalues[:] = np.array([0] + list(hist_fuck.binvalues[1:]))


hist_n.Project(tree_n, "neutrons.hits.tof")
hist_ph.Project(tree_ph, "photons.hits.tof", max_events=1E5)
hist_n *= 1.0/hist_n.bin_width/n_pulses_n

hist_ph += hist_fuck




hist_ph *= hist_n.binvalues[hist_n.FindBin(30)] /hist_ph.binvalues[-1]
hist_n += np.concatenate([hist_ph.binvalues, np.zeros(len(hist_n) - len(hist_ph))])



hist_n.Draw()
hist_ph.Draw()
mt2.thesis_plot([hist_n])




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
