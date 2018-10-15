import ROOT
import numpy as np
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F
import os


treeSP, n_pulsesSP = mt2.NCorrRun('SP', "Cf252").neutrons_doubles_tree

nbins = 65

histSP = TH1F(0,10,nbinss=nbins)
hist_theory = TH1F(0,10,nbinss=nbins)

# histErg = TH1F(0,10,nbinss=nbins*10)

E = np.array(hist_theory.bincenters[0])
# histErg.binvalues = np.e**(-0.88*E)*np.sinh( (2.0*E)**0.5)
erg = np.e**(-0.88*E)*np.sinh( (2.0*E)**0.5)
# histErg.__update_hist_from_containers__()

# for i in range(600000):
#     hist_theory.Fill(histErg.GetRandom())
# hist_theory.update_bin_containers_from_hist()

cut = "neutrons.hits.ForwardDet == 0"
histSP.Project(treeSP, "neutrons.hits.erg",cut)

erg *= sum(histSP.binvalues)/sum(erg)

histSP /= erg
histSP.Draw()
hist_theory.Draw("same hist")
# histErg.Draw()



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






