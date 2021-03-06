
import numpy as np

import numpy as np

import pickle
import matplotlib as mpl
import sys
import matplotlib.cm as cm


from mpl_toolkits.mplot3d import Axes3D


mpl.use('TkAgg')
import mytools2 as mt2
import mytools as mt
import ROOT

ROOT.TGaxis.SetMaxDigits(2)
from TH1Wrapper import *
import mytools2 as mt2
#
# view = ROOT.TView.CreateView(1)
# view.SetRange(5,5,5,25,25,25);
# import ctypes
# irep = ctypes.c_int()
# view.SetView(90,0,90,irep);

forward = True

tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP= mt2.NCorrRun('DP','DU', Forward=True).neutrons_doubles_tree

# theta_cuts = np.linspace(24,180,6)
theta_cuts = [24,165,180]
theta_cuts = list(zip(theta_cuts[0:-1],theta_cuts[1:]))

def transpose(hist):
    _hist = hist.__copy__()
    hist.transpose()
    hist += _hist
    hist.Multiply(0.5)
    for i in range(len(hist.binvalues)):
        hist.binvalues[i][i+1:] = 0
    return hist

histos= []
n_erg_bins = 2
max_erg = 6


for i, thetas in enumerate(theta_cuts):


    cut = mt.cut_rangeAND(thetas, '180/3.1415*neutrons.coinc_hits[0].coinc_theta')
    if not forward:
        cut += '&& (neutrons.coinc_hits[0].ForwardDet[0]==0 && neutrons.coinc_hits[0].ForwardDet[1]==0)'
    histSP = TH2F(0, max_erg, n_erg_bins)
    histDP = TH2F(0, max_erg, n_erg_bins)
    drw = 'neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]'

    print ('SP, {0}: {1} events'.format(thetas, histSP.Project(tree_SP, drw, cut, weight=1.0/n_pulses_SP)))
    print ('DP, {0}: {1} events\n'.format(thetas, histDP.Project(tree_DP, drw, cut, weight=1.0/n_pulses_DP)))

    transpose(histSP)
    transpose(histDP)

    histSP -= 0.5*histDP
    histSP /= histDP
    histos.append(histSP)

bar_args = []

lines = []
c1 = ROOT.TCanvas()
c1.Divide(len(histos))

cd_i = 1
_max = 0

data = OrderedDict()

for theta_cut, histSP in zip(theta_cuts, histos):
    __ = np.array(list(itertools.product(histSP.__binLeftEdges__[0], histSP.__binLeftEdges__[1]))).transpose()
    _x = __[0]
    _y = __[1]
    w = histSP.bincenters[0][1] - histSP.bincenters[0][0]

    bar_args.append((_x, _y, np.zeros_like(_x), w, w, histSP.binvalues.flatten()))

    theta_cut = list(map(int,theta_cut))
    histSP.SetTitle('{0}^{{#circ}} #leq #theta_{{nn}}<{1}^{{#circ}}'.format(*theta_cut))

    _hist = histSP.__copy__()
    _hist *= 4

    c1.cd(cd_i)
    cd_i += 1
    ROOT.gPad.SetPhi(-30)

    histSP.Draw('Lego', make_new_canvas=False)
    #
    # histSP.SetTitleSize(.3);

    histSP.SetStats(0)

    histSP.GetXaxis().SetTitle('E2 [MeV]')
    histSP.GetXaxis().SetLabelSize(.06);

    histSP.GetYaxis().SetTitle('E1 [MeV]')
    histSP.GetYaxis().SetLabelSize(.06);

    histSP.GetZaxis().SetTitleSize(.05);
    histSP.GetZaxis().CenterTitle()
    histSP.GetZaxis().SetTitle('Y_{corr}/Y_{uncorr}')
    histSP.GetZaxis().SetLabelSize(.06);

    histSP.GetZaxis().SetTitleOffset(1.)

    histSP.GetYaxis().SetTitleOffset(1.5)
    histSP.GetXaxis().SetTitleOffset(1.5)

    if max(histSP.binvalues.flatten())>_max:
        _max = max(histSP.binvalues.flatten())
    print('theta cut', theta_cut)
    theta_cut = tuple(theta_cut)

    data[theta_cut] = []

    def fix(n):
        return float('{:.2f}'.format(n))
    for i in range(len(histSP.binvalues)):
        for j in range(0,i+1):
            line = ROOT.TPolyLine3D()
            x, y = histSP.bincenters[0][i], histSP.bincenters[1][j]
            z = histSP.binvalues[i][j]
            err = histSP.binerrors[i][j]

            data[theta_cut].append(
                ((fix(x),fix(y)), (fix(z), fix(err)))
            )

            line.SetPoint(0,x,y,z)
            line.SetPoint(1,x,y,z +err)
            line.Draw('same')
            print('bin {0:.0f},{1:.0f}; value = {3:.2f} +/- {2:.2f}'.format(x,y,err, z))
            lines.append(line)
print(data)
for hist in histos:
    hist.SetMaximum(1.15*_max)





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
