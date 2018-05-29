import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2
import ROOT
import numpy as np
import mytools2 as mt2
import mytools as mt
import pickle

target = "Al"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_singles_tree

photon_dict = []
for degree in mt.angles:
    hist = TH1F(-35,150,binwidths=2,title=degree)
    hist.Project(treeSP_doubles, 'neutrons.hits.tof', 'neutrons.hits.det == {}'.format(degree))
    photon_dict.append(hist.binvalues)
    hist.Draw()

file = open("Al.pickle", 'w')
pickle.dump(photon_dict, file)


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
