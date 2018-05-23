import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2

import ROOT
import numpy as np
import mytools2 as mt2

target = "DU"

# treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
#

hist = TH1F(9-10,9+15,binwidths=0.5)

for i in np.random.poisson(1, 20000):
    hist.Fill(2*i + np.random.randn(1) + 7.5)


hist.GetXaxis().SetTitle("PMT timing sum")
hist.GetYaxis().SetTitle("counts")
hist.Draw('E')
ROOT.gStyle.SetOptStat('rm')
hist.SetTitle('')
hist.SetLineWidth(2)
mt2.thesis_plot(hist)

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
