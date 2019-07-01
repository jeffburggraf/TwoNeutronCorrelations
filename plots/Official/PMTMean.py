import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2

import ROOT
import numpy as np
import mytools2 as mt2

from scipy import special

target = "DU"

# treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
#

hist = TH1F(-2, 10,binwidths=0.5)
hist_p = TH1F(-2, 10,binwidths=0.15)

mu = 2.6

for i, x in enumerate(hist_p.bincenters[0]):
    y = np.e**(-mu)*mu**x/special.factorial(x, exact=False)
    hist_p.binvalues[i] = y

hist_p.__update_hist_from_containers__()

for i in range(10000):
    hist.Fill(hist_p.GetRandom() + np.random.randn()*0.9)

hist /= hist.bin_width
hist.GetXaxis().SetTitle("mean PMT #Deltat [ns]")
hist.GetYaxis().SetTitle("total counts [s^{-1}]")
hist.Draw('E')
ROOT.gStyle.SetOptStat('rm')
hist.SetTitle('')
hist.SetLineWidth(2)
ROOT.TGaxis.SetMaxDigits(3)

print(hist)
mt2.thesis_plot(hist, True)

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
