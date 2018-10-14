import sys
sys.path.append('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis')
import ROOT
from TH1Wrapper import TH1F
import numpy as np
import mytools2 as mt2
import mytools as mt


np.random.seed(3)
import matplotlib as mpl
font = {'family':'DejaVu Sans',
        'size': 17}
mpl.use('TkAgg')
mpl.rc('font', **font)

mpl.rc('text', usetex=True)

import matplotlib.pyplot as plt

# ====================
binwidth = 12
n_erg_bins = 4
plot_FREYA = True
DP_max_events = 600000
correction = True
avg_cut = True
# ====================

fig = None

# erg_bins = [0,1.5,3,4.5,9]

__tree__, _ = mt2.NCorrRun('SP', "Cf252").neutrons_doubles_tree
__erg_hist__ = TH1F(0.4, 9, 100)
__erg_hist__.Project(__tree__, "neutrons.coinc_hits.erg")

erg_bins, _ = mt2.median(__erg_hist__, n_erg_bins)
cut_off_bins = [0.4,1.4,2.4,3.4]

markers = ["<","^","o","v","s","D"]

assert len(markers) >= n_erg_bins

def gen_plots(target, plot, smooth = 1):
    FREYA = "FREYA" in target

    global binwidth, n_erg_bins, fig, c1, erg_bins, DP_max_events

    rebin_factor = 4
    _min_bin = 0

    nbins =int((180-24)/binwidth)

    if smooth:
        nbins *= rebin_factor

    print("getting trees")
    treeSP, n_pulsesSP = mt2.NCorrRun('SP',target).neutrons_doubles_tree
    treeDP, n_pulsesDP = mt2.NCorrRun('DP',target).neutrons_doubles_tree
    print("got trees")

    if DP_max_events is not None:
        treeDP.GetEntry(DP_max_events)
        _n_pulsesDP = treeDP.PulseNumber

        if _n_pulsesDP == 0:
            DP_max_events = None
        else:
            n_pulsesDP = _n_pulsesDP

    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(6, 11.5))
    histos_new = []
    _max_new = 0

    drw = "180/3.1415*neutrons.coinc_hits[0].coinc_theta"

    if correction:
        cut_DP = get_weighted_cut(target, subtract_accidentals=False, max_events=1000000)

    # for index, (E1, E2) in enumerate(zip(erg_bins[0:-1], erg_bins[1:])):
    for index, cut_value in enumerate(cut_off_bins if not avg_cut else zip(erg_bins[0:-1], erg_bins[1:])):
        E1=E2=0

        if avg_cut:
            E1, E2 = cut_value
            cut = mt.cut_rangeAND([E1, E2], "0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])")
        else:
            E_cut = cut_value
            cut = "(neutrons.coinc_hits[0].erg[0]>{0} && neutrons.coinc_hits[0].erg[1]>{0})".format(E_cut)

        histSP = TH1F(_min_bin, 180, nbinss=nbins)
        histDP = TH1F(_min_bin, 180, nbinss=nbins)
        histDP_weighted = TH1F(_min_bin, 180, nbinss=nbins)

        histDP_weighted.Project(treeDP, drw, cut="{0}*({1})".format(cut, cut_DP), max_events=DP_max_events if not FREYA else None)
        histDP_weighted /= n_pulsesDP
        histSP.Project(treeSP, drw, cut=cut, weight=1.0/n_pulsesSP)
        print("project 1 done for {}".format(target))
        histDP.Project(treeDP, drw, cut=cut, weight=1.0/n_pulsesDP, max_events=DP_max_events if not FREYA else None)
        print("project 2 done for {}\n".format(target))

        histDP_weighted *= sum(histDP.binvalues)/sum(histDP_weighted.binvalues)

        if smooth:
            histSP = histSP.MySmooth(smooth*rebin_factor, int(rebin_factor))
            histDP = histDP.MySmooth(smooth*rebin_factor, int(rebin_factor))
            histDP_weighted = histDP_weighted.MySmooth(smooth * rebin_factor, int(rebin_factor))

        histSP.SetMinimum(0)

        if correction:
            histSP /= (0.5*histDP_weighted)
        else:
            histSP /= (0.5*histDP)

        if not avg_cut:
            title_mpl = r" E>{0:.1f}$".format(E_cut)
        else:
            title_mpl = r"${0:.1f}<\overline{{ E_{{n}} }}<{1:.1f}$".format(E1, E2)

        if plot:
            # Note: Because this histogram is really a Kernal Density Estimte, the error bars are really a "fill between" region)
            ax.errorbar(histSP.bincenters[0], histSP.binvalues, yerr=histSP.binerrors, linewidth=1
                        , elinewidth=1., mec='black', capsize=2, c='black', drawstyle='steps-mid', label=title_mpl,
                        marker=markers[index])  # drawstyle='steps-mid'
            ax.set_xlabel(r"$\theta_{nn}$")
            if index == 0:
                ax.set_ylabel(r"$nn_{corr}/nn_{uncorr}$")

        # test for highest value
        if max(histSP.binvalues)>_max_new:
            _max_new = max(histSP.binvalues)

        histos_new.append(histSP)

    if plot:
        _max = 0
        for hist in histos_new:
            if max(hist.binvalues)>_max:
                _max = max(hist.binvalues)

        ax.grid()
        ax.minorticks_on()

        _ticks = map(lambda x: float("{:.1E}".format(x)), np.linspace(0,1.2*_max,5))
        ax.set_yticks(_ticks)
        ax.set_ylim(0, 1.2 * _max)
        ax.set_xlim(0, 190)
        ax.set_xticks(np.arange(0, 180 + 30, 30))

        plt.subplots_adjust(top=0.98, bottom=.08, wspace=0.07, hspace=0.30, left=0.15, right=0.97)
        plt.minorticks_on()

        ROOT.gSystem.ProcessEvents()

    return histos_new, (ax if plot else None)

DU_histos_l, ax  = gen_plots("Cf252", 1)
Freya_histos_l, _ = gen_plots("FREYA_Cf252", 0, 0.5)

if plot_FREYA:
    index = 0
    for i, (histDU, histF) in enumerate(zip(DU_histos_l, Freya_histos_l)):
        ax.errorbar(histF.bincenters[0], histF.binvalues * sum(histDU.binvalues)/sum(histF.binvalues),
            linestyle ='--', marker=markers[index], fillstyle='none')
        index += 1

plt.legend(handlelength=0)
ax.legend(handlelength=0)

plt.savefig("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/{}.png"
            .format("FinalCF252Result_w_FREYA({})".format("thresh_cut" if not avg_cut else "avg_cut")))
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

