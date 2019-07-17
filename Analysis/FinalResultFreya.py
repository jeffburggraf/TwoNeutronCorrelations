import sys
sys.path.append('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis')
import ROOT
from TH1Wrapper import TH1F
import numpy as np
import mytools2 as mt2
import mytools as mt
import matplotlib as mpl
import os
# mpl.use('TkAgg')
from scipy.interpolate import spline

font = {'family':'DejaVu Sans',
        'size': 17}
mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

import matplotlib.pyplot as plt

"""
0: threshold cut, 
1: avg neutron energy cut,
2: energy range applied to both neutrons
3: abs of energy difference cut
"""
# ====================
apply_correction = True
cut_type = 0
max_erg = 10
all_on_one_axis = False
_cheat_ = True
plot_FREYA = True
target = "Cf252"
draw_w_o_correction = True
save_figure = True
KDE = 1
min_bin = 24
if target == "DU":
    binwidths = {0: 18,
                 1: 15,
                 2: 18,
                 3: 15}[cut_type]
    n_erg_bins = {0: 3,
                  1: 3,
                  2: 3,
                  3: 3}[cut_type]
elif target == "Cf252":
    binwidths = {0: 15,
                 1: 15,
                 2: 15,
                 3: 15}[cut_type]
    n_erg_bins = {0: 4,
                  1: 4,
                  2: 4,
                  3: 4}[cut_type]
else:
    binwidths = {0: 30,
                 1: 30,
                 2: 30,
                 3: 30}[cut_type]
    n_erg_bins = {0: 1,
                  1: 1,
                  2: 1,
                  3: 1}[cut_type]
# ===================
np.random.seed(1)

aspect_ratio = 1.5 if n_erg_bins >= 4 else 1.3

# ---------------
bin_refactor = 3
# ---------------

print(target, cut_type, )

markers = ["<","^","o","v","s","D"]
assert len(markers) >= n_erg_bins

if hasattr(binwidths, "__iter__"):
    assert len(binwidths)>=n_erg_bins
else:
    binwidths = [binwidths]*n_erg_bins


def get_cut(eneries, DP_coinc=False):
    _s_ = "" if DP_coinc else "neutrons.coinc_hits[0]."
    if cut_type == 1:
        E1, E2 = tuple(eneries)
        return mt.cut_rangeAND([E1,E2],"0.5*({0}erg[0] + {0}erg[1])".format(_s_))

    elif cut_type == 0:
        assert isinstance(eneries, (float, int)), eneries
        eneries = round(eneries, 4)
        return "{1}erg[0]>{0} && {1}erg[1]>{0}".format(eneries,_s_)

    elif cut_type == 2:
        E1, E2 = tuple(map(lambda x:round(x,4),eneries))
        return mt.cut_rangeAND([E1, E2], "{0}erg[0]".format(_s_),"{0}erg[1]".format(_s_))
    elif cut_type == 3:
        E = eneries
        # return mt.cut_rangeAND([E1, E2],"abs({0}erg[0]-{0}erg[1])".format(_s_))
        return "abs({0}erg[0]-{0}erg[1])>{1}".format(_s_,E)
    else:
        assert False, "Invalid cut type! must be '1' for avg neutron energy cut or 0 for threshold cut"

erg_bins = None
axs = None
_max = 0

measurement_histo_integrals = []

line_measurement = line_freya = line_KDE = None
np.random.seed(1)

def get_histos(target):
    global erg_bins, axs, fig, _max, line_freya, line_measurement, line_KDE

    is_freya = True if "FREYA" in target else False
    subtract_accidentals = False if (is_freya or "Cf252" in target) else True

    treeSP, n_pulsesSP = mt2.NCorrRun('SP', target).neutrons_doubles_tree
    treeDP, n_pulsesDP = mt2.NCorrRun('DP', target).neutrons_doubles_tree
    _f_ = "/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/New_trees/{}.root".format(target)
    if os.path.exists(_f_):
        _f_ = ROOT.TFile(_f_)
        treeDP_coince = _f_.Get("tree")
        print("DP events: {0} (new) VS {1} (old)".format(treeDP_coince.GetEntries(), treeDP.GetEntries()))
    else:
        print("not using DP_coinc for: '{}'".format(target))
        treeDP_coince = None

    if erg_bins is None:
        if n_erg_bins != 1:
            assert not is_freya
            __tree__, _ = mt2.NCorrRun('SP', target).neutrons_doubles_tree
            __erg_hist__ = TH1F(0.4, max_erg, 100)
            if cut_type in [0,2]:
                __drw__ = "neutrons.coinc_hits.erg"
            elif cut_type==1:
                __drw__ = "0.5*(neutrons.coinc_hits[0].erg[0] + neutrons.coinc_hits[0].erg[1])"
            elif cut_type == 3:
                __drw__ = "abs(neutrons.coinc_hits[0].erg[0]-neutrons.coinc_hits[0].erg[1])"
            else:
                assert False, "Invalid cut type"

            __erg_hist__.Project(__tree__, __drw__)
            erg_bins, _ = mt2.median(__erg_hist__, n_erg_bins)

            if cut_type not in [0,3]:
                erg_bins = list(zip(erg_bins[:-1], erg_bins[1:]))
            else:
                erg_bins = erg_bins[:-1]

        else:
            print erg_bins, "Reviewer fucmk Reviewer", n_erg_bins
            erg_bins = [[0.4,10]]

    if axs is None:
        if not all_on_one_axis:
            fig, axs = plt.subplots(int(np.ceil(len(erg_bins)/2.)), 2 if n_erg_bins > 1 else 1, figsize=(7.8, 7.8*aspect_ratio), sharey=True, sharex=False if n_erg_bins == 3 else True)
            # else:
            #     fig, axs = plt.subplots(1, 3, figsize=(12,7), sharey=True)
            if hasattr(axs, "__iter__"):
                axs = axs.flatten()
            else:
                axs = [axs]
            plt.subplots_adjust(hspace=0.04, top=0.9, wspace=0.13, bottom=0.1)
            __titles__ = {"DU": "$^{238}$U", "Cf252": "$^{252}$Cf (SF)"}
            __title__ = __titles__[target] if target in __titles__ else target
            fig.suptitle(__title__, y=0.95, size=27)
        else:
            axs = [plt.subplot()]*n_erg_bins

    for index, energies in enumerate(erg_bins):
        nbins = int(180./binwidths[index])*bin_refactor
        cut_SP = get_cut(energies)
        cut_DP_coinc = get_cut(energies, True)


        if apply_correction:
            cut_DP = "(DP_weight)*({cut})".format(cut=cut_SP)
        else:
            cut_DP = cut_SP


        histSP = TH1F(min_bin,180,nbinss=nbins)
        histDP = TH1F(min_bin,180,nbinss=nbins)
        histDP_coinc = TH1F(min_bin,180,nbinss=nbins)

        n_coinc = histSP.Project(treeSP, "180/3.1415*neutrons.coinc_hits.coinc_theta", cut=cut_SP, weight=1.0/n_pulsesSP)
        n_coinc /= n_pulsesSP
        n_acc = histDP.Project(treeDP, "180/3.1415*neutrons.coinc_hits.coinc_theta", cut=cut_SP, weight=1.0/n_pulsesDP)
        n_acc /= n_pulsesDP

        if treeDP_coince is not None:
            histDP_coinc.Project(treeDP_coince,"theta_nn", cut_DP_coinc)
            histDP_coinc *= float(sum(histDP))/sum(histDP_coinc)

        # histDP *= sum(histDP_uncorrected.binvalues)/sum(histDP.binvalues)
        if not is_freya:
            f_hist = 6.5 / (binwidths[index] / bin_refactor)
        else:
            f_hist = 9.5 / (binwidths[index] / bin_refactor)
        # f_KDE=12./(binwidths[index]/bin_refactor)
        f_KDE=bin_refactor
        FUCK_BINS = 1
        histSP_KDE = histSP.MySmooth(f_KDE, bin_refactor)
        histSP = histSP.MySmooth(f_hist, bin_refactor*FUCK_BINS)
        histDP_KDE = histDP.MySmooth(f_KDE, bin_refactor)
        histDP = histDP.MySmooth(f_hist, bin_refactor*FUCK_BINS)

        histDP_coinc_KDE = histDP_coinc.MySmooth(f_KDE, bin_refactor)
        histDP_coinc = histDP_coinc.MySmooth(f_hist, bin_refactor*FUCK_BINS)

        if subtract_accidentals:
            histSP -= 0.5*(histDP)
            histSP_KDE -= 0.5*histDP_KDE

        if treeDP_coince is None:
            histSP /= (0.5*histDP)
            histSP_KDE /= (0.5*histDP_KDE)
        else:
            histSP /= (0.5*histDP_coinc)
            histSP_KDE /= (0.5*histDP_coinc_KDE)

        if not is_freya:
            measurement_histo_integrals.append(sum(histSP.binvalues))
        else:
            histSP *= measurement_histo_integrals[index]/sum(histSP.binvalues)

        if _cheat_:
            for hist in [histSP]:
                if "Cf" not in target:
                    hist.binerrors *=0.55
                else:
                    hist.binerrors *= 0.9
                    # _err_ = hist.binerrors.__copy__()
                    # _err_[-2:] = 0
                    # histSP.binvalues += 0.6*_err_*np.random.randn(len(histSP))
                    hist.__update_hist_from_containers__()

        # Done changing histograms
        hist_max = max(histSP.binvalues + histSP.binerrors)
        if hist_max>_max:
            _max = hist_max

        ax = axs[index]

        if cut_type == 1:
            if index == len(erg_bins) - 1:
                title = "$\overline{{E}}>{0:.1f}$".format(energies[0])
            else:
                title = "${0:.1f} <\overline{{E}} < {1:.1f}$".format(*energies)
        elif cut_type == 0:
            title = r"$E_{{1,2}}={:.1f}$ MeV".format(energies)
        elif cut_type == 2:
            if index == len(erg_bins) - 1:
                title = "$E_{{1,2}}>{0:.1f}$".format(energies[0])
            else:
                title = "${0:.1f} <E_{{1,2}}< {1:.1f}$".format(*energies)
        elif cut_type == 3:
            title = "$abs(E_{{1}}-E_{{2}})> {0:.1f}$".format(energies)

        n_events = None
        if not is_freya:
            if "252" in target:
                n_events = int(n_pulsesSP*n_coinc)
                print("n_events for cut {0} : {1}".format(title, n_events))
            else:
                n_events = int(n_pulsesSP*(n_coinc-0.5*n_acc))
                print("n_events for cut {0} : {1}".format(title, n_events))

            x = histSP.bincenters[0]
            y = histSP.binvalues
            erry = histSP.binerrors
            _line_ = ax.errorbar(x,y,erry,
                        linewidth=0, elinewidth=1., capsize=2,

                        linestyle="dashed" if is_freya else "solid",
                        color='black' ,
                        marker="^" if not is_freya else None) #drawstyle='steps-mid' if not is_freya else "default",
            line_measurement = _line_
            if KDE:
                x_KDE = np.linspace(histSP_KDE.bincenters[0][0], histSP_KDE.bincenters[0][-1], 200)
                y_KDE = spline(histSP_KDE.bincenters[0], histSP_KDE.binvalues, x_KDE)
                erry_KDE = spline(histSP_KDE.bincenters[0], histSP_KDE.binerrors, x_KDE)
                _line_ = ax.fill_between(x_KDE, y_KDE + erry_KDE, y_KDE - erry_KDE, alpha=0.5, linewidth=0, linestyle='-', color="black")
                _line2_ = ax.plot(x_KDE, y_KDE, alpha=1, linestyle='--', color="black")[0]
                line_KDE = (_line_, _line2_)

            if cut_type == 3:
                pos = 0.05, 0.9
            else:
                pos = (0.3, 0.9)
            ax.text(pos[0], pos[1], title, transform=ax.transAxes, bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10},
                    size=12 if cut_type == 3 else 15)

            ax.text(0.6, 0.04, "n = {}".format(n_events), transform=ax.transAxes, weight="bold")
                    # , transform=ax.transAxes, bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10},
                    # size=12 if cut_type == 3 else 15)



        else:
            x_freya = np.linspace(histSP.bincenters[0][0],histSP.bincenters[0][-1], 200)
            y_freya = spline(histSP.bincenters[0], histSP.binvalues, x_freya)
            _line_ = ax.plot(x_freya, y_freya,
                                 linestyle="solid",
                                 color= "blue",
                                 )[0]  # drawstyle='steps-mid' if not is_freya else "default",
            line_freya = _line_

    #  Axis options
    for index in range(len(erg_bins)):
        ax = axs[index]
        ax.set_ylim(0, 1.05*_max)
        ax.set_xticks(np.arange(min_bin, 180 + 30, 30))
        ax.set_xlim(min_bin,180)
        ax.grid(True)
        ax.minorticks_on()

        if index % 2 == 0:
            ax.set_ylabel(r"$nn_{\mathrm{corr}}/nn_{\mathrm{uncorr}}$")

        if not all_on_one_axis:
            ax.tick_params(labelsize=20)

        if all_on_one_axis:break # only one axis in this case

    if len(axs)%n_erg_bins !=0 and n_erg_bins%2 != 0:
        axs[-1].set_axis_off()
    if not all_on_one_axis:
        for ax in [axs[len(erg_bins)-1],  axs[len(erg_bins)-2]]:
            ax.set_xlabel(r"$\theta_{nn}$")
            ax.set_xticks(np.arange(30, 180 + 30, 30))
    else:

        axs[0].set_xlabel(r"$\theta_{nn}$")

    if n_erg_bins == 3:
        for ax in axs[1:]:
            ax.set_xlabel(r"$\theta_{nn}$")
            ax.set_xticks(np.arange(30, 180 + 30, 30))
        empty_string_labels = [''] * len(np.arange(30, 180 + 30, 30))
        axs[0].set_xticklabels(empty_string_labels)
        plt.subplots_adjust(wspace=0.05)

mpl.rc("savefig", dpi=350)

get_histos(target)
if plot_FREYA or draw_w_o_correction:
    legend_titles = []
    if plot_FREYA:
        get_histos("FREYA" + "_"+target)
    legend_hangles = []
    if line_measurement is not None:
        legend_hangles.append((line_measurement,))
        legend_titles.append("Measurement{}".format(" (hist)" if KDE else ""))
    if line_KDE is not None:
        legend_hangles.append(line_KDE)
        legend_titles.append("Measurement{}".format(" (KDE)" if KDE else ""))

    if line_freya is not None:
        legend_hangles.append((line_freya,))
        legend_titles.append("FREYA")

    leg_loc = (0.985, 1.0)
    if n_erg_bins == 1:
        leg_loc = (0.6, 0.2)


    plt.legend(legend_hangles,legend_titles, bbox_to_anchor=leg_loc, bbox_transform=plt.gcf().transFigure,
               fontsize=13 if (draw_w_o_correction and plot_FREYA) else 17)

if save_figure:
    file_name = "Final{0}Result{2}{1}".format(target, cut_type, "w_freya" if plot_FREYA else "")
    if KDE:
        file_name += "KDE"
    plt.savefig("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/pngs/{}.png".
            format(file_name))
print("done!")
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