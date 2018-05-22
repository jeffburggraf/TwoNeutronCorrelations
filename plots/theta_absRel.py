import ROOT
import numpy as np
import mytools2 as mt2
from TH1Wrapper import TH1F
import mytools as mt

targets = "DU", 'Pu'

# subtract = 'U233'

treeSP1, pulses_SP1_doubles = mt2.NCorrRun("SP", targets[0], generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeSP2, pulses_SP2_doubles = mt2.NCorrRun("SP", targets[-1], generate_dictionary=False,Forward=True).neutrons_doubles_tree


# erg_bins =mt2.median(D2OErgHist, 4)[0]
# if erg_bins[0]<0.5:
#     erg_bins[0] = 0.5

# dx_pos = (30*2.5)/3.
# binszpos = np.arange(-2.5*15, 2.5*15 +dx_pos, dx_pos , dtype=np.float32)

binsphi = []
for degree in mt.angles:
    binsphi.append(degree - 10)
    binsphi.append(degree + 10)
    if degree>180:break

binsphi = [0] + binsphi

hist1 = TH1F(binarrays=binsphi)
hist2 = TH1F(binarrays=binsphi)

hist1.Project(treeSP1, '180/3.1415*neutrons.coinc_hits.theta_abs')
hist2.Project(treeSP2, '180/3.1415*neutrons.coinc_hits.theta_abs')


hist1 /= hist2

hist1.SetTitle('{0}/{1}'.format(*targets))
hist1.GetXaxis().SetTitle("#theta_{abs}")
hist1.GetYaxis().SetTitle("ratio")

hist1.Draw()
mt2.thesis_plot(hist1)


if __name__ == "__main__":
    import ROOT as ROOT
    from multiprocessing import Process, Queue
    import time, sys, os


    def input_thread(q, stdin):
        while True:
            print 'ROOT: '
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
                        print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()


