import ROOT
import mytools2 as mt2
import numpy as np
from matplotlib import pyplot as plt
from TH1Wrapper import TH1F
import mytools as mt
from scipy.special import erfc
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#  ===================
target = "DU"
draw_type = "divide"
bin_width = {"lowest": 0.75, "highest": 0.75, "divide":0.15, "asym":0.15}[draw_type]
min_erg = {"lowest": 0.35, "highest": 0.35, "divide":0, "asym":0}[draw_type]
max_erg = {"lowest": 5, "highest": 8, "divide":1, "asym":0.9}[draw_type]
theta_bins = [100, 150, 170, 180]
plot_p_value = True
#  ===================

theta_bins = list(zip(theta_bins[:-1], theta_bins[1:]))  # + [(30, 180)]

draw_type = draw_type.lower()
tree_Acc, n_pulses_Acc = mt2.NCorrRun("DP","DU",).neutrons_doubles_tree
tree_SP, n_pulsesSP = mt2.NCorrRun("SP","DU",).neutrons_doubles_tree

if target == "DU":
    _f = ROOT.TFile("/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/new_trees/DU.root")
    tree_DP = _f.Get("tree")
else:
    assert target == "Cf252"
    tree_DP = tree_Acc


def get_draw(tree, theta_cut):
    if hasattr(tree, "neutrons"):
        E1 = "neutrons.coinc_hits[].erg[0]"
        E2 = "neutrons.coinc_hits[].erg[1]"
        theta = "neutrons.coinc_hits[].coinc_theta*180/3.14"
        print tree

    else:
        E1 = "erg[0]"
        E2 = "erg[1]"
        theta = "theta_nn"
        print tree

    if draw_type == "lowest":
        result = "{E1}*({E1}<{E2}) + {E2}*({E2}<{E1})".format(E1=E1,
                                                            E2=E2)
    elif draw_type == "highest":
        result = "{E1}*({E1}>{E2}) + {E2}*({E2}>{E1})".format(E1=E1,
                                                            E2=E2)
    elif draw_type == "divide":
        result = "{E1}/{E2}*({E1}<{E2}) + {E2}/{E1}*({E2}<{E1})".format(E1=E1,
                                                              E2=E2)
    elif draw_type == "asym":
        result = "abs({E1}-{E2})/({E1} + {E2})".format(E1=E1, E2=E2)
    else:
        assert False

    cut = mt.cut_rangeAND(theta_cut, theta)
    print result, cut
    return result, cut

markers = ["o", "v", "d","x","P"]


ys = []
yerrs = []
for index, thetas in enumerate(theta_bins):
    hist_SP = TH1F(min_erg, max_erg, binwidths=bin_width)
    hist_DP = TH1F(min_erg, max_erg, binwidths=bin_width)
    hist_Acc = TH1F(min_erg, max_erg, binwidths=bin_width)

    hist_SP.Project(tree_SP, *get_draw(tree_SP, thetas), weight=1.0/n_pulsesSP)
    good_bins = np.where(hist_SP.binvalues>0)
    hist_DP.Project(tree_DP, *get_draw(tree_DP, thetas))
    hist_Acc.Project(tree_Acc, *get_draw(tree_Acc, thetas), weight=1.0/n_pulses_Acc)

    hist_Acc *= 0.5

    hist_DP *= sum(hist_Acc.binvalues)/sum(hist_DP.binvalues)

    hist_SP -= hist_Acc

    hist_SP /= hist_DP

    x = np.array(hist_SP.bincenters[0])[good_bins]
    y = hist_SP.binvalues[good_bins]
    yerr = hist_SP.binerrors[good_bins]

    ys.append(y)
    yerrs.append(yerr)

    plt.errorbar(x, y, yerr, label=r"${0}^{{\circ}}<\theta_{{nn}}<{1}^{{\circ}}$".format(*thetas), marker=markers[index])


    # hist_SP.Draw()





if draw_type == "highest":
    x_label = "Energy of the most energetic neutron in pair [MeV]"
elif draw_type == "lowest":
    x_label = "Energy of the least energetic neutron in pair [MeV]"
elif draw_type == "divide":
    x_label = "Ratio of neutron with smaller energy to neutron with larger energy"
elif draw_type == "asym":
    x_label = "Normalized neutron energy difference"

y_label = "$nn_{corr}/nn_{uncorr}$"

plt.ylabel(y_label, fontsize=15)
plt.xlabel(x_label, fontsize=15)
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc="lower left")
plt.subplots_adjust(bottom=0.15)


plt.show()


# hist_SP.GetXaxis().SetTitle(x_label)
# hist_SP.GetYaxis().SetTitle(y_label)
# mt2.thesis_plot([hist_SP])
# hist_SP.SetTitle("#theta_{{nn}} > {0} degrees".format(theta_cut_off))
# hist_SP.SetStats(0)

tb = ROOT.TBrowser()

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









