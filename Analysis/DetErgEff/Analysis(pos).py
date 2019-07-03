import ROOT
import numpy as np
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F
import os
from matplotlib import pyplot as plt

nbins= 5
max_events = None
n_plots = 2,2

treeSP, n_pulsesSP = mt2.NCorrRun('SP', "Cf252").neutrons_doubles_tree


class Det:
    def __init__(self, angle, top=None):
        if top is not None:
            assert angle in [30, 330]
            self.cut = "neutrons.hits.det == {} && ForwardDet == {topbot}".format(angle, topbot=1 if top is True else -1)
            self.label = "{0} ({top})".format(angle, top="top" if top is True else "bottom")
            self.hist = TH1F(-40,40, nbinss=1)
            self.c = 10
            self.x = [19.05] if top is True else [-19.05]
            self.xerr = 2

        else:
            self.cut = "neutrons.hits.det == {}".format(angle)
            self.label = "{0}".format(angle)
            self.hist = TH1F(-40, 40, nbinss=nbins)
            self.c = self.hist.bin_width
            self.x = self.hist.bincenters[0]
            self.xerr = 0.7*self.hist.bin_width


dets = []
for angle in mt.angles:

    if angle in [30, 330]:pass
        # dets.append(Det(angle, True))
        # dets.append(Det(angle, False))
    else:
        dets.append(Det(angle))


histos = []

means = []

_min = None
_max = None
for det in dets:
    hist = det.hist
    hist.Project(treeSP, "neutrons.hits.z", det.cut, max_events=max_events if max_events is not None else treeSP.GetEntries()-1)
    hist.update_bin_containers_from_hist()
    # hist = hist.MySmooth(1, rebin=1)
    hist /= det.c

    means.append(np.mean(hist.binvalues))
    histos.append(hist)

    if _max is None:
        _max = max(hist.binvalues)
        _min = min(hist.binvalues)
    else:
        _min = min([_min, min(hist.binvalues)])
        _max = max([_max, max(hist.binvalues)])

hist_mean = np.mean(means)

fig, axs = plt.subplots(*n_plots, sharey=True, sharex=True, figsize=(12,6))
axs = axs.flatten()

_n_plots = n_plots
n_plots = n_plots[0]*n_plots[1]


line_styles = ["-","--","-."]

_min *= 1.0 / hist_mean
_max *= 1.0 / hist_mean

for index, det in enumerate(dets):
    det.hist /= hist_mean
    # det.xerr /= hist_mean
    hist = det.hist
    print hist.binvalues, det.label

    ax_i = index % (n_plots)
    ax = axs[ax_i]
    line_i = index//n_plots
    ax.errorbar(det.x, hist.binvalues, yerr=det.hist.binerrors, marker="^", markersize=3, xerr=det.xerr, elinewidth=0.75, linewidth=1.2, linestyle=line_styles[line_i%len(line_styles)], label=det.label)

    leg = ax.legend(title="det. angle", loc="upper right", handlelength=3)
    # for leg in leg.legendHandles:
    #     leg.set_linewidth(5.0)

    ax.set_ylim(0.6*_min, 1.1*_max)
    ax.grid(True)

    # if ax_i % (_n_plots[0]) == 0:
    #     ax.set_label("rel. neutron efficiency ")

y_title = "rel. neutron efficiency "
x_title = "vertical position [cm]"

fig.text(0.07, 0.5, y_title , va='center', rotation='vertical')
fig.text(0.45, 0.05, x_title, va='center')

fig_ = plt.figure()


confidence_half_with = np.std([det.hist.binvalues for det in dets], axis=0)
eff_mean = np.mean([det.hist.binvalues for det in dets], axis=0)
eff_mean_err = np.sqrt(np.mean([det.hist.binerrors**2 for det in dets], axis=0))


plt.fill_between(hist.bincenters[0], eff_mean + confidence_half_with, eff_mean - confidence_half_with, alpha=0.5, color="grey", linewidth=0, label=r"$\pm$ the S.D. of all detectors")
plt.errorbar(hist.bincenters[0], eff_mean, yerr=eff_mean_err, linestyle="--", color="black", linewidth=0.6, marker="^")

plt.legend(loc="lower center")
plt.ylim(0)

plt.xlabel(x_title)
plt.ylabel(y_title)

plt.subplots_adjust(wspace=0.03, hspace=0.03)

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






