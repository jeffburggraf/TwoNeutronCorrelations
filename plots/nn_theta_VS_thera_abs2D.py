
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

from TH1Wrapper import *
import mytools2 as mt2


forward = True

tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP= mt2.NCorrRun('DP','DU', Forward=True).neutrons_doubles_tree

# theta_cuts = [24,100,150,180]
theta_cuts = [24,120,160,180]
theta_cuts = list(zip(theta_cuts[0:-1],theta_cuts[1:]))


histos= []
theta_abs_bins = []
for angle in mt.angles:
    theta_abs_bins.extend([angle-12, angle+12])
    if angle == 150:
        break

for i, thetas in enumerate(theta_cuts):
    cut = mt.cut_rangeAND(thetas, '180/3.1415*neutrons.coinc_hits[0].coinc_theta')
    if not forward:
        cut += '&& (neutrons.coinc_hits[0].ForwardDet[0]==0 && neutrons.coinc_hits[0].ForwardDet[1]==0)'
    histSP = TH2F(binarrays=theta_abs_bins)
    histDP = TH2F(binarrays=theta_abs_bins)

    # histSP = TH1F(binarrays=theta_abs_bins)
    # histDP = TH1F(binarrays=theta_abs_bins)
    drw = '180/3.1415*neutrons.coinc_hits[0].theta_abs[0]:180/3.1415*neutrons.coinc_hits[0].theta_abs[1]'
    # drw = '180/3.1415*neutrons.coinc_hits[0].theta_abs[1]'

    print ('SP, {0}: {1} events'.format(thetas, histSP.Project(tree_SP, drw, cut, weight=1.0/n_pulses_SP)))
    print ('DP, {0}: {1} events\n'.format(thetas, histDP.Project(tree_DP, drw, cut, weight=1.0/n_pulses_DP)))

    histSP = 0.5*(histSP + histSP.transpose())
    histDP = 0.5*(histDP + histDP.transpose())

    histSP -= 0.5*histDP
    histSP /= histDP
    histos.append(histSP)

    histSP.set_max_rel_err(0.3)

bar_args = []

lines = []
c1 = ROOT.TCanvas()
c1.Divide(len(histos))

cd_i = 1
_max = 0

pads = []

data = OrderedDict()

for theta_cut, histSP in zip(theta_cuts, histos):
    theta_cut = list(map(int,theta_cut))
    histSP.SetTitle('{0}^{{#circ}} #leq #theta_{{nn}}<{1}^{{#circ}}'.format(*theta_cut))

    _hist = histSP.__copy__()
    _hist *= 4

    pad = c1.cd(cd_i)
    pads.append(pad)
    cd_i += 1
    ROOT.gPad.SetPhi(-30)

    histSP.Draw('Lego', make_new_canvas=False)
    histSP.SetMinimum(0)


    # histSP.GetYaxis().SetNoExponent(ROOT.kTRUE);
    # histSP.GetXaxis().SetNoExponent(ROOT.kTRUE);

    if max(histSP.binvalues.flatten())>_max:
        _max = max(histSP.binvalues.flatten())

    theta_cut = tuple(theta_cut)

    data[theta_cut] = []

    def fix(n):
        return float('{:.2f}'.format(n))

    for i in range(len(histSP.binvalues)):

        for j in range(0, i+1):
            line = ROOT.TPolyLine3D()
            x, y = histSP.bincenters[0][i], histSP.bincenters[1][j]
            z = histSP.binvalues[i][j]
            err = histSP.binerrors[i][j]

            data[theta_cut].append( ((fix(x),fix(y)), (fix(z), fix(err))) )

            line.SetPoint(0, x, y, z)
            line.SetPoint(1, x, y, z+err)

            line.Draw('same')
            lines.append(line)

axis_label_size = 0.09

for hist,pad in zip(histos,pads):
    hist.GetXaxis().SetNdivisions(4, 5, 0, 0);
    hist.GetYaxis().SetNdivisions(4, 5, 0, 0);

    pad.cd()
    hist.SetMaximum(1.05*_max)
    ROOT.gStyle.SetTitleFontSize(0.1)


    hist.GetZaxis().SetTitle('nn_{corr}/nn_{uncorr}')
    hist.GetZaxis().CenterTitle()
    hist.GetZaxis().SetLabelSize(axis_label_size)
    hist.GetZaxis().SetTitleSize(.08)
    hist.GetZaxis().SetTitleOffset(1.35)

    hist.GetXaxis().SetTitle('#theta1_{abs}')
    hist.GetXaxis().SetTitleSize(.1)
    hist.GetXaxis().SetTitleOffset(1.)
    hist.GetXaxis().CenterTitle()

    hist.GetXaxis().SetLabelOffset(.035)
    hist.GetXaxis().SetLabelSize(axis_label_size)

    hist.GetYaxis().SetTitle('')
    hist.GetYaxis().SetTitleSize(.1)
    hist.GetYaxis().SetTitleOffset(1.1)
    hist.GetYaxis().CenterTitle()

    hist.GetYaxis().SetLabelSize(axis_label_size)

    pad.SetLeftMargin(0.2)
    pad.SetRightMargin(0.16)
    pad.SetBottomMargin(0.2)

    hist.SetStats(0)

    # hist.SetTitleSize(0.2)

hist = TH1F(0,1,10)
hist.GetXaxis().SetTitle('#theta2_{abs}')
hist.Draw()

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
