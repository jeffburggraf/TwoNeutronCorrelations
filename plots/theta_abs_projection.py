import numpy as np
from TH1Wrapper import TH1F,TH2F
import mytools2 as mt2
import os
from itertools import product
import wpca
import ROOT
import mytools as mt
#
tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP = mt2.NCorrRun('DP','DU', Forward=True).neutrons_doubles_tree
#

theta_abs_bins = []
for angle in mt.angles:
    theta_abs_bins.extend([angle-10, angle+10])
    if angle == 150:
        break

#  ===============
theta_range = [150,180]
# =============

histSP = TH2F(binarrays=theta_abs_bins)
histDP = TH2F(binarrays=theta_abs_bins)

cut =mt.cut_rangeAND(theta_range, 'neutrons.coinc_hits[].coinc_theta*180/3.1415')

histSP.Project(tree_SP, '180/3.141*neutrons.coinc_hits[0].theta_abs[0]:180/3.141*neutrons.coinc_hits[0].theta_abs[1]', cut=cut, weight=1.0/n_pulses_SP)
histDP.Project(tree_DP, '180/3.141*neutrons.coinc_hits[0].theta_abs[0]:180/3.141*neutrons.coinc_hits[0].theta_abs[1]', cut=cut, weight=1.0/n_pulses_DP)

c1 = ROOT.TCanvas('c1','c1',700, 700)

histSP -= 0.5* histSP

histSP /= histDP

histSP = 0.5*(histSP.transpose(copy=True) + histSP)

histSP.Draw('colz', make_new_canvas=0)
histSP.set_max_rel_err(0.5)


# bin_indicies = np.array(list(product(range(len(_hist.bincenters[0])), range(len(_hist.bincenters[1])))))
n = len(histSP.bincenters[0])
bin_indicies = (list(product(range(n), range(n))))
X = []
weights = []



for I, (i1, i2) in enumerate(bin_indicies):
    center = [histSP.bincenters[0][i1], histSP.bincenters[1][i2]]
    weight = np.sqrt(abs(histSP.binvalues[i1,i2]))

    if weight>0:
        X.append(center)
        weights.append([weight, weight])


X = np.array(X)
weights = np.array(weights)


pca = wpca.WPCA()
pca.fit_transform(X, weights=weights)

mean = np.mean(X, axis = 0)

ars = []
for i in range(2):
    c = pca.components_[i]
    s = np.sqrt(pca.explained_variance_[i])
    x1, y1 = tuple(mean)
    x2, y2 = tuple(mean + s*c)
    print(x1, y1, x2, y2)
    a = ROOT.TArrow(x1, y1, x2, y2, 0.03)
    if i == 0:
        a.SetLineColor(ROOT.kRed)
    ars.append(a)
    a.Draw()


proj_hist = histSP.OneDProjection(mean, pca.components_[1], 10, True)
proj_hist.MySmooth(1)
proj_hist.Draw()
print(np.sqrt(pca.explained_variance_))


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



