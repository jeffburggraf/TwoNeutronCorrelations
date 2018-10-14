import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt

from scipy.stats import gaussian_kde

import  matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

binwidth = 12
# f = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/DP_FREYA_neutrons_coinc.root ")
# tree = f.Get("tree")

treeSP, n_pulsesSP = mt2.NCorrRun('SP',"DU").neutrons_doubles_tree
treeDP, n_pulsesDP = mt2.NCorrRun('DP',"DU").neutrons_doubles_tree

histSP = TH1F(24,180,binwidths=0.05)
histDP = TH1F(24,180,binwidths=0.05)

drw = "180/3.1415*neutrons.coinc_hits.coinc_theta"
histSP.Project(treeSP, drw,  weight =1.0/n_pulsesSP )
histDP.Project(treeDP, drw, weight=1.0/n_pulsesDP)


# histSP.Draw()
pointsSP = (np.array(histSP.bincenters[0])[np.where(histSP.binvalues!=0)])
pointsDP = (np.array(histDP.bincenters[0])[np.where(histDP.binvalues!=0)])


class KDE_result:
    def __init__(self, x,y,erry):
        self.x = x
        self.y = y
        self.erry = erry

    def __div__(self, other):
        if isinstance(other, KDE_result):
            assert len(self) == len(other)
            erry = np.sqrt((other.erry ** 2 * self.y** 2 + self.erry ** 2 * other.y** 2) / (other.y**4))
            return KDE_result(self.x, self.y/other.y,erry)

    def __len__(self):
        return len(self.x)

class KDE(KDE_result):
    def __init__(self, hist, sigma):
        assert isinstance(hist, TH1F)
        width = int(len(hist)*(2/3.))
        _x_ = np.linspace(-width,width, 3*len(hist))
        kernel = np.e**(-_x_**2/(2.*sigma**2))
        kernel = kernel/sum(kernel)

        data = np.concatenate((hist.binvalues[::-1], hist.binvalues, hist.binvalues[::-1]))
        y = np.convolve(kernel, data,mode="same")[len(hist):-len(hist)]
        erry = np.sqrt(y)
        x = hist.bincenters[0]

        KDE_result.__init__(self,x,y,erry)


# kde_SP = gaussian_kde(pointsSP)
# kde_DP = gaussian_kde(pointsDP)

kde_SP = KDE(histSP, 20)
kde_DP = KDE(histDP, 20)
# kde.fit(points)

r = kde_SP/kde_DP

X = np.linspace(0,180,200)
plt.plot(r.x, r.y)
# plt.plot(kde_DP.x, kde_DP.y)
plt.show()




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




