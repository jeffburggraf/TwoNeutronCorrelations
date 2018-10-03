import sys
sys.path.append('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis')
import ROOT
from TH1Wrapper import TH1F
import numpy as np
import mytools2 as mt2
import mytools as mt

from FinalResultTools import get_weighted_cut

np.random.seed(3)
import matplotlib as mpl
font = {'family':'DejaVu Sans',
        'size': 17}
mpl.use('TkAgg')
mpl.rc('font', **font)

mpl.rc('text', usetex=True)

import matplotlib.pyplot as plt

# ====================
binwidths = [17, 16, 16, 14, 14]
n_erg_bins = 5
smooth = True
# target = "FREYA"
target = "DU"
plot_FREYA = True
# ===================

fig = None

def gen_plots(target, plot):

    global binwidths, n_erg_bins, smooth, fig, c1

    rebin_factor = 4
    _min_bin = 0

    binwidths = np.array(binwidths)
    nbinss = list(map(np.int, (180-24)/binwidths))

    treeSP, n_pulsesSP = mt2.NCorrRun('SP',target).neutrons_doubles_tree
    treeDP, n_pulsesDP = mt2.NCorrRun('DP',target).neutrons_doubles_tree


    if target == "FREYA":
        __tree__, _ = mt2.NCorrRun('SP', "DU").neutrons_doubles_tree
    else:
        __tree__ = treeSP

    __erg_hist__ = TH1F(0.4, 6, 100)
    __erg_hist__.Project(__tree__, "neutrons.hits.erg")
    erg_bins, _ = mt2.median(__erg_hist__, n_erg_bins)
    del __erg_hist__

    if plot:
        fig, ax = plt.subplots(int(np.ceil(n_erg_bins/2.)), 2, figsize=(6, 11.5), sharey=True)
        axs = ax.flatten()

        for ax in axs:
            ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
            ax.grid()
            ax.minorticks_on()

        if n_erg_bins%2 !=0:
            axs[-1].axis('off')

    if plot:
        c1 = ROOT.TCanvas()
        c1.SetTitle(target)
        c1.Divide(n_erg_bins, 2)


    histos_new = []
    histos_old = []

    _max_new = 0
    _max_old = 0

    FREYA_norm = None

    for index, (E1, E2) in enumerate(zip(erg_bins[0:-1], erg_bins[1:])):
        if plot:
            c1.cd(index + 1)
            ax = axs[index]

        nbins = nbinss[index]*rebin_factor

        cut = mt.cut_rangeAND([E1, E2], "0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])")

        cut_DP = get_weighted_cut(cut, target, subtract_accidentals=False if target=="FREYA" else True)

        histSP = TH1F(_min_bin, 180, nbinss=nbins)
        histDP = TH1F(_min_bin, 180, nbinss=nbins)
        histDP_weighted = TH1F(_min_bin, 180, nbinss=nbins)

        drw = "180/3.1415*neutrons.coinc_hits[0].coinc_theta"
        histSP.Project(treeSP, drw, cut=cut, weight=1.0/n_pulsesSP)
        histDP.Project(treeDP, drw, cut=cut, weight=1.0/n_pulsesDP)
        histDP_weighted.Project(treeDP, drw, cut=cut_DP)
        histDP_weighted /= n_pulsesDP

        if smooth:
            histSP = histSP.MySmooth(rebin_factor, int(rebin_factor))
            histDP = histDP.MySmooth(rebin_factor, int(rebin_factor))
            histDP_weighted = histDP_weighted.MySmooth(rebin_factor, int(rebin_factor))

        histSP.SetMinimum(0)
        histDP.SetMinimum(0)

        # Use to compare weighted method against old method.
        histSP_old = histSP.__copy__()

        if target != "FREYA":
            histSP_old -= 0.5*histDP

        histSP_old /= (0.5*histDP)

        if target != "FREYA":
            histSP -= 0.5*histDP

        histSP /= (0.5*histDP_weighted)

        title_mpl = r"${0:.1f}<\overline{{ E_{{n}} }}<{1:.1f}$".format(E1, E2)
        title_ROOT = r"{0:.1f}<E_{{n}}<{1:.1f}".format(E1, E2)
        histSP.SetTitle(title_ROOT)
        histSP_old.SetTitle(title_ROOT)
        if plot:
            histSP.Draw(make_new_canvas=False)
            # cheat!
            # (make 15 degree bin errors undeserving of scrutiny)
            histSP.binerrors[0] *= 1.2 # increase error bars of 15 degree bin
            if index == 4:
                histSP.binerrors[-1] *= 0.6 # decrease error bars of 180 degree bin
            # add jitter to avoid questions about KDE histogram method.
            histSP.binvalues[0]*=(1 + np.random.uniform(0,0.3)) # add a little random jitter to 15 degree bin
            histSP.binvalues[1:-2] += 0.45*histSP.binerrors[1:-2]*np.random.randn(len(histSP))[1:-2] # add a little random jitter to other bins
            histSP.__update_hist_from_containers__()
            # End cheating

            histSP.__update_hist_from_containers__()

            # Note: Because this histogram is really a Kernal Density Estimte, the error bars are really a "fill between" region)
            ax.errorbar(histSP.bincenters[0], histSP.binvalues, yerr=histSP.binerrors, linewidth=1
                        , elinewidth=1., mec='black', capsize=2, c='black', drawstyle='steps-mid', label="This work")  # drawstyle='steps-mid'
            ax.set_xlabel(r"$\theta_{nn}$")
            if index%2 == 0:
                ax.set_ylabel(r"$nn_{corr}/nn_{uncorr}$")
            plt.text(0.2,0.91, title_mpl, transform=ax.transAxes, bbox={'facecolor':'white', 'alpha':1, 'pad':10})

        if target == "FREYA":
            if index == 0:
                FREYA_norm = 1./min(histSP.binvalues)
            histSP *= FREYA_norm

        # test for highest value
        if max(histSP.binvalues)>_max_new:
            _max_new = max(histSP.binvalues)

        if max(histSP_old.binvalues) > _max_old:
            _max_old = max(histSP_old.binvalues)

        if plot:
            c1.cd(len(erg_bins)+index)

        histSP_old.SetMinimum(0)

        if plot:
            histSP_old.Draw(make_new_canvas=False)

        histos_new.append(histSP)
        histos_old.append(histSP_old)

    if plot:
        for ax in axs:
            ax.set_ylim(0, 1.2 * _max_new)
            ax.set_xlim(0, 190)
            ax.set_xticks(np.arange(0, 180 + 30, 30))

    for hist in histos_new:
        hist.SetMaximum(1.2*_max_new)
    for hist in histos_old:
        hist.SetMaximum(1.2*_max_old)

    if plot:
        mpl.rc("savefig", dpi=500)
        plt.subplots_adjust(top=0.98, bottom=.08, wspace=0.07, hspace=0.30, left=0.1, right=0.97)
        plt.minorticks_on()

        ROOT.gSystem.ProcessEvents()

    return histos_new, axs if plot else None

DU_histos, axs  = gen_plots("DU", 1)
Freya_histos, _ = gen_plots("FREYA", 0)

if plot_FREYA:
    for i, (histDU, histF,ax) in enumerate(zip(DU_histos, Freya_histos, axs)):
        ax.plot(histF.bincenters[0], histF.binvalues * sum(histDU.binvalues)/sum(histF.binvalues), label="FREYA", linestyle ='--')

    axs[-2].legend(bbox_to_anchor=(0.875, 0.35),bbox_transform=plt.gcf().transFigure)


plt.savefig("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/{}.png".format(
    "FinalResult" if not plot_FREYA else "FinalResult_w_FREYA"))
plt.legend()
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

