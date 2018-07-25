import ROOT
import numpy as np
import mytools as mt
from collections import OrderedDict

import mytools2 as mt2

from TH1Wrapper import TH1F

from scipy.optimize import curve_fit
from numpy.polynomial.legendre import legfit
import pickle


forward = True


treeSP, n_pulses_SP = mt2.NCorrRun("SP", 'DU', generate_dictionary=False, Forward=True).neutrons_doubles_tree
treeDP, n_pulses_DP = mt2.NCorrRun("DP", 'DU', generate_dictionary=False, Forward=True).neutrons_doubles_tree

_cut_ = '(neutrons.coinc_hits[0].ForwardDet[0]==0 && neutrons.coinc_hits[0].ForwardDet[1] == 0)'

min_erg = 0.4
erg_hist_SP = TH1F(min_erg,5.5,binwidths=0.001)
erg_hist_DP = TH1F(min_erg,5.5,binwidths=0.001)
erg_hist_SP.Project(treeSP, '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])', '' if forward else _cut_, weight=1.0/n_pulses_SP)
erg_hist_DP.Project(treeDP, '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])', '' if forward else _cut_, weight=1.0/n_pulses_DP)

raw_energy = erg_hist_SP.__copy__()
erg_hist_SP -= 0.5*erg_hist_DP

erg_hist_DP.binerrors = raw_energy.binerrors
erg_hist_DP *= sum(raw_energy.binvalues)/sum(erg_hist_DP.binvalues)
erg_bins, _ = mt2.median(erg_hist_SP,4)
# w = 1.
# erg_bins = np.arange(0.5,5*w, w)



NN = 300
erg_hist_SP = erg_hist_SP.Rebin(NN)
raw_energy = raw_energy.Rebin(NN)
erg_hist_DP = erg_hist_DP.Rebin(NN)

raw_energy.SetLineColor(ROOT.kRed)


erg_hist_SP.Draw()
raw_energy.Draw('same')

histos = []

c1 = ROOT.TCanvas()
c1.Divide(len(erg_bins)-1,1)

global_max = 0

def Legendre(cos_X, *p):
    return np.l

data_dict = {}
data_dict['cuts'] = OrderedDict()

data_dict['energy'] = {'SP':{'x':erg_hist_SP.x, 'y':erg_hist_SP.binvalues, 'err':erg_hist_SP.binerrors},
                       'DP': {'x': erg_hist_DP.x, 'y': erg_hist_DP.binvalues, 'err': erg_hist_DP.binerrors},
                       'raw': {'x': raw_energy.x, 'y': raw_energy.binvalues, 'err': raw_energy.binerrors}}

#################################
bin_width = 15
courseness = float(5)  # Arg to hist.Rebin
Smooth = [0.6,0.73, 0.68, 0.5]
method = 'pad'
#################################

data_dict['theta_bin_width'] = np.mean(Smooth)*bin_width


for i, (e0, e1) in enumerate(zip(erg_bins[0:-1], erg_bins[1:])):
    if isinstance(Smooth, list):
        smooth = Smooth[i]
    else:
        smooth = Smooth
    c1.cd(i+1)
    histSP = TH1F(24, 180, binwidths=bin_width/courseness)
    histDP = TH1F(24, 180, binwidths=bin_width/courseness)
    erg_cut = '0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])'
    erg_cut = mt.cut_rangeAND([e0, e1], erg_cut) + '&& neutrons.coinc_hits.ForwardTopBot == 0'

    if not forward:
        erg_cut += ' && ' + _cut_

    drw = '180/3.1415*neutrons.coinc_hits.coinc_theta'
    n_coinc = histSP.Project(treeSP, drw, erg_cut, weight=1.0/n_pulses_SP)
    n_accidentals = histDP.Project(treeDP,drw, erg_cut, weight=1.0/n_pulses_DP)
    n_accidentals = int(0.5*n_accidentals/n_pulses_DP*n_pulses_SP)

    if smooth:
        histSP.MySmooth(courseness*smooth/2, edge_method=method)
        histDP.MySmooth(courseness*smooth/2, edge_method=method)

    I0 = np.sum(histSP.binvalues)
    histSP -= 0.5*histDP
    I1 = np.sum(histSP.binvalues)

    n_events = int(n_pulses_SP*sum(histSP.binvalues))

    histSP = histSP.Rebin(courseness)
    histDP = histDP.Rebin(courseness)

    histSP /= histDP
    if smooth:
        histSP.MySmooth(smooth,edge_method=method)

    histSP.SetMinimum(0)

    print('Between {0:.2f} and {1:.2f}:\nn_trues: {2}'.format(e0,e1, n_events))
    print('n_accidentals: {0}'.format(n_accidentals))
    print('n_coinc: {0}'.format(n_coinc))
    print('Accidental subtraction change: {:.2f}'.format((I1-I0)/I0))


    if i == 0:
        norm_factor = 1./max(histSP.binvalues)

    histSP *= norm_factor
    title = '{0:.1f}<E<{1:.1f}'.format(e0,e1)
    histSP.SetTitle(title)

    if max(histSP) > global_max:
        global_max = max(histSP)

    histos.append(histSP)

    hist_x_data = np.array([histSP.GetBinCenter(i)for i in range(1,len(histSP.binvalues) + 1)])
    hist_y_data = np.array([histSP.GetBinContent(i)for i in range(1,len(histSP.binvalues) + 1)])
    hist_y_derr = np.array([histSP.GetBinError(i)for i in range(1,len(histSP.binvalues) + 1)])
    fit_x_data = np.cos(hist_x_data*3.1415/180)

    best_params, diagnostic = legfit(fit_x_data, hist_y_data, 3, w=1.0 / hist_y_derr, full=True)
    print(diagnostic)
    print('Dipole-monopole ratio: {0:.2f}\n'.format(best_params[2]/best_params[0]))
    data_dict['x_data'] = hist_x_data

    data_dict['cuts'][title] = {'params': best_params}
    data_dict['cuts'][title]['param_errors'] = diagnostic
    data_dict['cuts'][title]['y_data'] = hist_y_data
    data_dict['cuts'][title]['y_err'] = hist_y_derr
    data_dict['cuts'][title]['cut_range'] = [e0,e1]
    data_dict['cuts'][title]['n_trues'] = n_events
    data_dict['cuts'][title]['n_accidentals'] = n_accidentals
    data_dict['cuts'][title]['n_coinc'] = n_coinc

f = open('FinalDUResultData.pickle','wb')
pickle.dump(data_dict, f)
f.close()

##################   ROOT  #################_
norm_factor = sum(histos[0].binvalues)
for i, hist in enumerate(histos):
    c1.cd(i+1)
    hist.SetMaximum(global_max*1.1)

    hist.SetMinimum(0)
    hist.GetXaxis().SetRange(0, 180);
    hist.Draw(make_new_canvas=False)


# tb = ROOT.TBrowser()


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
#
