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


#   ============
erg_bin_width = 0.5
theta_bins = np.linspace(120,180, 4)[:-1]
cut_type = "max"
subtract_accidentals = True
two_event_uncorr = True
#    =========
max_erg = 5 if cut_type=="diff" else 3
min_erg = 0.75 if cut_type == "max" else 0.5

tree_SP, n_pulses_SP = mt2.NCorrRun('SP','DU', Forward=True).neutrons_doubles_tree
tree_DP, n_pulses_DP = mt2.NCorrRun('DP','DU', Forward=True).neutrons_doubles_tree

fig, ax = plt.subplots(1,1)
# axs = axs.flatten()

cut_type = cut_type.lower()

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

def get_n_cut(energy, angle, two_events):
    assert cut_type in ["min", "max", "diff"]
    if two_events is False:
        erg_str = "neutrons.coinc_hits[].erg"
        theta_str = "neutrons.coinc_hits[].coinc_theta*180/3.14"
    else:
        erg_str = "erg"
        theta_str = "theta_nn"

    if cut_type in ["min", "max"]:
        result =  "{erg_str}[0]{gtlt}{0} && {erg_str}[1]{gtlt}{0}".\
            format(energy, gtlt=">" if cut_type == "min" else "<", erg_str=erg_str)
    else:
        result = "abs({erg_str}[0] - {erg_str}[1])>{0}".format(energy, erg_str=erg_str)

    return "{result} && {theta_str}>{angle}".format(result=result, theta_str=theta_str, angle=angle)


markers = ["o", "v", "d","x","P"]

_f = ROOT.TFile("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/new_trees/DU.root")
tree_uncorr_two_events = _f.Get("tree")

for index, theta in enumerate(theta_bins):
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

        cut = get_n_cut(e, theta, False)
        cut_two_events = get_n_cut(e, theta, True)
        # print cut
        # print cut_two_events
        # assert False

        n_SP = tree_SP.GetEntries(cut)
        n_DP = tree_DP.GetEntries(cut)
        n_DP_two_events = tree_uncorr_two_events.GetEntries(cut_two_events)

        y_SP.append(n_SP)
        y_DP.append(n_DP)
        y_DP_two_events.append(n_DP_two_events)

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

    label = r"$\theta_{{nn}} > {0}^{{\circ}}$".format(theta)
    # label = r"{0}$^{{\circ}}$< $\theta_{{nn}}$ < {1}$^{{\circ}}$".format(*theta_range)

    ax.errorbar(x, y.y, yerr=y.yerr, label=label, capsize=5, marker=marker, linestyle="--", linewidth=0.8)
    ax.set_ylabel(y_label, fontsize=14)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xticks(x_ticks)
    x_label = {"max":"max. neutron energy", "min":"min. neutron energy", "diff":"min. absolute neutron energy difference"}[cut_type]
    ax.set_xlabel(x_label + " [MeV]", fontsize=13)
    # ax.text(0.06, 0.16, "n$_{{2}}$ energy threshold = {0} MeV".format(fixed_erg_cut_max),transform=ax.transAxes)
    # ax.text(0.06, 0.1, r"{0}$^{{\circ}}$< $\theta_{{nn}}$ < {1}$^{{\circ}}$".format(min_angle, max_angle), transform=ax.transAxes)
    ax.legend(fontsize=13)
    ax.grid()


plt.show()



