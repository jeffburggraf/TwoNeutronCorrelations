import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2

hist =TH1F(0,18,binwidths=1)

for i in np.random.randn(20000):
    i = i*3.2 + 9
    hist.Fill(i)

hist -= 30

hist.Draw("")
hist.SetName("")
hist.SetMinimum(0)
hist.SetMarkerStyle(32)

ROOT.gStyle.SetOptStat("mr")

hist.SetStats(1)

hist.GetXaxis().SetTitle("PMT timing average")
hist.GetYaxis().SetTitle("counts")

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
