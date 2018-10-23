import ROOT
import numpy as np
from TH1Wrapper import TH1F
import mytools as mt
import mytools2 as mt2

hist = TH1F(4,10.5,binwidths = 0.1)
hist_test = TH1F(4,10.5,binwidths = 0.1)

xs = mt.ENDF_crossection("/Users/jeffreyburggraf/FREYA/data_freya/92238(g,f)xs.txt")


for i in range(len(hist)):
    x = hist.bincenters[0][i]
    hist.binvalues[i] = np.e**(-0.54)*xs.Eval(x)


hist.__update_hist_from_containers__()

hist.normalize()

hist.GetXaxis().SetTitle("Energy of fission inducing photon [MeV]")
hist.GetYaxis().SetTitle("rel. rate [arb. units]")

hist.Draw()
hist.SetLineWidth(2)
hist.SetLineColor(ROOT.kBlack)

mt2.thesis_plot([hist], big_font=0.07)

ROOT.gStyle.SetOptStat('m')

path = "/Users/jeffreyburggraf/FREYA/MyFREYA/2nCorr/G_ergs.txt"
f = open(path, "w")

for i in range(1000000):
    if i:
        f.write("\n")
    f.write("{:.2f}".format(hist.GetRandom()))
f.close()

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
