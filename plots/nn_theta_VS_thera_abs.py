
import numpy as np

import numpy as np

import pickle
import matplotlib as mpl
import sys
import matplotlib.cm as cm


from mpl_toolkits.mplot3d import Axes3D


forward  = False
mpl.use('TkAgg')
import mytools2 as mt2
import mytools as mt
import ROOT

# ROOT.TGaxis.SetMaxDigits(3)
from TH1Wrapper import *
import mytools2 as mt2

case = None
############################
extenuate_zero = False
extenuate_90 = False

range_90 = [78, 105]
range_zero = [30, 60], [115, 180]

if extenuate_zero and not extenuate_90:
    n_bins = 12
    smooth = 0
    bin_factor = 3
    case = 1
    range_zero = [30, 55], [180-55, 180]
elif not extenuate_zero and extenuate_90:
    n_bins = 12
    smooth = 2
    bin_factor = 3
    case = 2
    range_90 = [79, 101]
else:
    n_bins = 14
    smooth = 2.
    bin_factor = 3
    case = 0


tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP= mt2.NCorrRun('DP','DU', Forward=True).neutrons_doubles_tree

theta_abs_hist = TH1F(0,180,binwidths=0.5)
theta_abs_hist.Project(tree_DP, '180/3.141*neutrons.coinc_hits.theta_abs')
theta_abs_hist.Draw('hist')


cut_90 = mt.cut_rangeOR(range_90, 'neutrons.coinc_hits[0].theta_abs[0]*180/3.14', 'neutrons.coinc_hits[0].theta_abs[1]*180/3.14')


    # not_cut_90 = '! (' + cut_90 + ')'

    # if case == 1:
__range_zero__ = [range_zero[0][1],range_zero[1][0]]
not_cut_90 = mt.cut_AND(mt.cut_range_outer(__range_zero__, 'neutrons.coinc_hits[0].theta_abs[0]*180/3.14'),
                        mt.cut_range_outer(__range_zero__, 'neutrons.coinc_hits[0].theta_abs[1]*180/3.14'))

if not forward:
    _ = ' && neutrons.coinc_hits[0].ForwardDet[0] == 0 && neutrons.coinc_hits[0].ForwardDet[1] == 0'
    cut_90 +=  _
    not_cut_90 += _

min_bin = 24
__n_bins__ = n_bins*bin_factor
if case == 1:
    min_bin = 100
    __n_bins__ = int(bin_factor * n_bins*0.5)
elif case == 2:
    pass

histSP_90 = TH1F(min_bin, 180, __n_bins__)
histDP_90 = TH1F(min_bin, 180, __n_bins__)


histSP_not_90 = TH1F(min_bin, 180, __n_bins__)
histDP_not_90 = TH1F(min_bin, 180, __n_bins__)

c1 = ROOT.TCanvas()
# c1.Divide(2)
# c1.cd(1)

histSP_90.Project(tree_SP, '180/3.1415*neutrons.coinc_hits[0].coinc_theta', weight=1.0/n_pulses_SP,cut=cut_90)
histDP_90.Project(tree_DP, '180/3.1415*neutrons.coinc_hits[0].coinc_theta', weight=1.0/n_pulses_DP, cut=cut_90)

histSP_90 -= 0.5*histDP_90

if smooth:
    histSP_90.MySmooth(smooth)
    histDP_90.MySmooth(smooth)

histSP_90 = histSP_90.Rebin(bin_factor)
histDP_90 = histDP_90.Rebin(bin_factor)

histSP_90 /= (0.5*histDP_90)
histSP_90.SetMarkerStyle(33)
histSP_90.SetMarkerSize(1.5)


histSP_90.Draw('hist E', make_new_canvas=0)
histSP_90.GetXaxis().SetTitle('#theta_{nn}')

histSP_90.GetYaxis().SetTitle('(nn_{corr})/(nn_{uncorr})')

var = '#theta_{abs}_{1,2}'
title_90 = (['#theta_{{abs}}_{{1}} *or* #theta_{{abs}}_{{2}} #in {rng}'.format(rng=range_90)]*2 + [''])[case]
# histSP_90.SetTitle(title_90)
histSP_90.SetLineStyle(7)
histSP_90.SetLineWidth(3)
mt2.thesis_plot([histSP_90], big_font=0.06)
histSP_90.GetXaxis().SetNdivisions(6,5,0,0);

c1.Modified()
c1.Update()
histSP_90.SetStats(0)

# c1.cd(2)

histSP_not_90.Project(tree_SP, '180/3.1415*neutrons.coinc_hits[0].coinc_theta', weight=1.0/n_pulses_SP,cut=not_cut_90)
histDP_not_90.Project(tree_DP, '180/3.1415*neutrons.coinc_hits[0].coinc_theta', weight=1.0/n_pulses_DP, cut=not_cut_90)

histSP_not_90 -= 0.5*histDP_not_90

if smooth:
    histSP_not_90.MySmooth(smooth)
    histDP_not_90.MySmooth(smooth)

histSP_not_90 = histSP_not_90.Rebin(bin_factor)
histDP_not_90 = histDP_not_90.Rebin(bin_factor)

histSP_not_90 /= (0.5*histDP_not_90)

_max = max(np.concatenate([histSP_not_90.binvalues, histSP_90.binvalues]))

histSP_not_90.Draw('hist E same',make_new_canvas=0)
histSP_not_90.GetXaxis().SetTitle('#theta_{nn}')
histSP_not_90.GetYaxis().SetTitle('(nn_{corr})/(nn_{uncorr})')
mt2.thesis_plot([histSP_not_90], big_font=0.05)
histSP_not_90.GetXaxis().SetNdivisions(6,5,0,0);
histSP_not_90.SetMarkerStyle(27)
histSP_not_90.SetLineWidth(3)
histSP_not_90.SetMarkerSize(1.5)


title_not_90= '#theta_{{abs}}_{{1}} *and* #theta_{{abs}}_{{2}} #in {0} #cup {1} '.format(*range_zero, var=var)
# histSP_not_90.SetTitle(title_not_90)
histSP_not_90.SetStats(0)


histSP_not_90.SetMinimum(0)
histSP_90.SetMinimum(0)
histSP_not_90.SetMaximum(round(1.5*_max, 1))
histSP_90.SetMaximum(round(1.4*_max, 1))

leg = ROOT.TLegend(0.2,0.75,0.88,0.95)
leg.AddEntry(histSP_90,     '{0}^{{#circ}}<#theta_{{abs}}<{1}^{{#circ}} for at least one neutron'.format(*range_90), 'lpf')
leg.AddEntry(histSP_not_90, '{0}^{{#circ}}<#theta_{{abs}}<{1}^{{#circ}} for neither neutrons'.format(*range_90), 'lpf')

leg.Draw()

c1.Modified()
c1.Update()

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
