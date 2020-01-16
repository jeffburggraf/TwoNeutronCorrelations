import numpy as np
from TH1Wrapper import TH1F,TH2F
from matplotlib import pyplot as plt
import mytools as mt
import mytools2 as mt2
import ROOT

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


#   ====================================
erg_bin_width = 0.5
theta_bins = [150,170, 170]
make_theta_bins_distinct = True
subtract_accidentals = True
two_event_uncorr = True
perp_cut_angle = 35  # None is not using
plot_p_value = False
save_fig = True
font_size = 10
#   ====================================

tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP = mt2.NCorrRun('DP','DU', Forward=True).neutrons_doubles_tree

# fig, ax = plt.subplots(1,1)
# axs = axs.flatten()

if make_theta_bins_distinct:
    if theta_bins[-1] != 180:
        theta_bins += [180]
    theta_bins = list(zip(theta_bins[:-1], theta_bins[1:]))


theta_bins = [theta_bins[0]] + [theta_bins[-1]]

def jitter(arr, p=0.04, max_n=4):
    if not hasattr(jitter, "n"):
        jitter.n = 0
    if jitter.n > max_n:
        jitter.n = 0

    n = jitter.n
    scale = (-1)**n*np.ceil(n/2.)
    jitter.n += 1
    print scale

    arr = np.array(arr)
    b_width = np.mean(arr[1:] - arr[:-1])
    # sigma = np.random.rand(len(arr))
    # sigma = np.where(abs(sigma)<3, sigma, 3)
    result = arr + p*b_width*scale
    return result

def get_n_cut(energy, angle, two_events, cut_type):
    assert cut_type in ["min", "max", "diff"]
    if two_events is False:
        erg_str = "neutrons.coinc_hits[].erg"
        theta_str = "neutrons.coinc_hits[].coinc_theta*180/3.14"
        theta_abs_str = "180/3.14*neutrons.coinc_hits[].theta_abs"
    else:
        erg_str = "erg"
        theta_str = "theta_nn"
        theta_abs_str = "180/3.14*abs_theta"

    if cut_type in ["min", "max"]:
        final_cut = "{erg_str}[0]{gtlt}{0} && {erg_str}[1]{gtlt}{0}".\
            format(energy, gtlt=">" if cut_type == "min" else "<", erg_str=erg_str)
    else:
        final_cut = "abs({erg_str}[0] - {erg_str}[1])>{0}".format(energy, erg_str=erg_str)

    if perp_cut_angle is not None:
        perp_cut = mt.cut_rangeOR([90-perp_cut_angle, 90+perp_cut_angle], "{0}[0]".format(theta_abs_str),"{0}[1]".format(theta_abs_str))
        final_cut = "{0} && ({1})".format(final_cut, perp_cut)

    if hasattr(angle,  "__iter__"):
        assert len(angle) == 2
        theta_nn_cut = mt.cut_rangeAND(angle, theta_str)
    else:
        theta_nn_cut = "{0}>{1}".format(theta_str, angle)

    final_cut = "{0} && {1}".format(final_cut, theta_nn_cut)

    return final_cut

markers = ["o", "v", "d","x","P"]

_f = ROOT.TFile("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/new_trees/DU.root")
tree_uncorr_two_events = _f.Get("tree")

fig, axs = plt.subplots(3,1, figsize=(5,8))
axs = axs.flatten()

for subplot_index, cut_type in enumerate(["diff", "max", "min"]):
    max_erg = 5.5 if cut_type == "diff" else 3
    min_erg = 0.75 if cut_type == "max" else 0.5
    ax = axs[subplot_index]

    XYs = []

    color_cycle = ["C{}".format(i) for i in range(10)]

    for index, theta in enumerate(theta_bins):
        if not make_theta_bins_distinct:
            theta = int(theta)
        # angle_cut = mt.cut_rangeAND(theta_range, "neutrons.coinc_hits[].coinc_theta*180/3.14")

        x = np.arange(min_erg, max_erg + erg_bin_width, erg_bin_width)
        y_SP = []
        y_DP = []
        y_DP_two_events = []

        marker = markers[index]

        for e in x:
            # other_n_cut0 = "neutrons.coinc_hits[].erg[0]<{}".format(e)
            # other_n_cut1 = "neutrons.coinc_hits[].erg[1]<{}".format(e)

            cut = get_n_cut(e, theta, False, cut_type)
            cut_two_events = get_n_cut(e, theta, True, cut_type)

            n_SP = tree_SP.GetEntries(cut)
            n_DP = tree_DP.GetEntries(cut)
            n_DP_two_events = tree_uncorr_two_events.GetEntries(cut_two_events)

            if n_SP>0:
                y_SP.append(n_SP)
                y_DP.append(n_DP)
                y_DP_two_events.append(n_DP_two_events)

        x = x[:len(y_SP)]
        if index == 0:
            print cut
            print cut_two_events

        y_SP = np.array(y_SP)
        y_DP = np.array(y_DP)
        y_DP_two_events = np.array(y_DP_two_events)

        y_SP = mt2.XYData(x, y_SP/n_pulses_SP, np.sqrt(y_SP)/n_pulses_SP)
        y_DP = mt2.XYData(x, y_DP/n_pulses_DP, np.sqrt(y_DP)/n_pulses_DP)
        y_DP_two_events = mt2.XYData(x, y_DP_two_events/n_pulses_DP, np.sqrt(y_DP_two_events)/n_pulses_DP)
        y_DP_two_events *= sum(y_DP.y)/sum(y_DP_two_events.y)

        if two_event_uncorr:
            y_DP = y_DP_two_events

        if subtract_accidentals is True:
            y_label = "nn$_{corr}$/nn$_{uncorr}$"
            y = (y_SP - y_DP/2.)/(y_DP/2.)
            if cut_type == "max":
                y += 0.3
        elif subtract_accidentals is False:
            y_label= "nn$_{SP}$/nn$_{uncorr}$"
            y = (y_SP)/(y_DP/2.)
        elif subtract_accidentals == "DP":
            y_label = "nn$_{DP}$"
            y = y_DP/2.
        elif subtract_accidentals == "SP":
            y_label = "nn$_{SP}$"
            y = y_SP
        else:
            assert False

        x_ticks = x[:]
        x = jitter(x)

        if not make_theta_bins_distinct:
            label = r"$\theta_{{nn}} > {0}^{{\circ}}$".format(theta)
        else:
            label = r"${0}^{{\circ}}<\theta_{{nn}} < {1}^{{\circ}}$".format(*theta)
        # label = r"{0}$^{{\circ}}$< $\theta_{{nn}}$ < {1}$^{{\circ}}$".format(*theta_range)

        ax.errorbar(x, y.y, yerr=y.yerr, label=label, capsize=4.5, marker=marker, markersize=4.5, linestyle="None", c=color_cycle[index])
        # ax.errorbar(x, y.y, yerr=0.4*y.yerr, label=label, capsize=4.5, marker='None', linestyle="None", c=color_cycle[index])
        ax.set_ylabel(y_label, fontsize=font_size*1.1)

        ax.tick_params(axis='x', labelsize=font_size*1.1)
        ax.tick_params(axis='y', labelsize=font_size*1.1)
        ax.set_xticks(x_ticks)
        x_label = {"max":"maximum energy of n-n pairs", "min":"minimum energy of n-n pairs", "diff":"minimum energy difference of n-n pairs"}[cut_type]
        ax.set_xlabel(x_label + " [MeV]", fontsize=font_size)
        # ax.text(0.06, 0.16, "n$_{{2}}$ energy threshold = {0} MeV".format(fixed_erg_cut_max),transform=ax.transAxes)
        # ax.text(0.06, 0.1, r"{0}$^{{\circ}}$< $\theta_{{nn}}$ < {1}$^{{\circ}}$".format(min_angle, max_angle), transform=ax.transAxes)
        ax.grid()

        if cut_type == "diff":
            ax.set_xticks(range(1,6))
        elif cut_type == "max":
            ax.set_xticks(np.arange(0.5,4, 0.5))
            ax.set_xlim(0.4, 3.75)


        XYs.append(y)

    if plot_p_value:
        ax2 = ax.twinx()

        color = "C3"

        ax2.plot(x + 0.05, mt2.p_value(XYs[0].y, XYs[0].yerr, XYs[-1].y, XYs[-1].yerr), linewidth=0, marker="o",
                 fillstyle='none', markersize=3.5, label="Compatibility test", color=color)

        # ax2.tick_params(axis='y', colors='black')
        # ax2.spines['right'].set_color('black')
        ax2.set_ylabel("p-value", rotation=90)
        ax2.set_yscale("log")
        ax2.set_yticks([1E0, 1E-1, 1E-2])

        ax2.spines['right'].set_color(color)
        ax2.yaxis.label.set_color(color)
        ax2.tick_params(axis='y', which='both', colors=color)
        ax2.legend(fontsize=font_size, bbox_to_anchor=(0.9, 0.92, 0, 0), bbox_transform=plt.gcf().transFigure, loc=4)

    # break

# ax2.legend(fontsize=font_size, bbox_to_anchor=(0.6, 3, 0.5,1))

# ax.legend(fontsize=font_size, bbox_to_anchor=(0, 3, 0.5,1), )
ax.legend(fontsize=font_size, bbox_to_anchor=(0.39, 0.92, 0,0), bbox_transform=plt.gcf().transFigure, loc=4)

plt.subplots_adjust(bottom=0.06, top=0.92, right=0.87, hspace=.35)

if save_fig:
    plt.savefig("/Users/jeffreyburggraf/Desktop/LargeAngleAnomaly.png", dpi=300)

plt.show()

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






