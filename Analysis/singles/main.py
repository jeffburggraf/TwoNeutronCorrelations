import ROOT
import numpy as np
from TH1Wrapper import  TH1F
import mytools2 as mt2
import mytools as mt

ROOT.gROOT.ProcessLine(".L twontree.h")

target = "Th"
max_events = None
double_only = False

treeCf, n_pulses_Cf = mt2.NCorrRun("SP", "Cf252").neutrons_doubles_tree
if double_only:
    tree_target, pulses_target = mt2.NCorrRun("SP", target).neutrons_doubles_tree
    tree_Al, pulses_Al = mt2.NCorrRun("SP", "Al").neutrons_doubles_tree
    tree_target_noise, pulses_target_noise = mt2.NCorrRun("SP", target).noise_doubles_tree
else:
    tree_target, pulses_target = mt2.NCorrRun("SP", target).neutrons_singles_tree
    tree_Al, pulses_Al = mt2.NCorrRun("SP", "Al").neutrons_singles_tree
    tree_target_noise, pulses_target_noise = mt2.NCorrRun("SP", target).noise_singles_tree

tree_Al_photons, pulses_Al_photon = mt2.NCorrRun("SP", "Al").photons_singles_tree
tree_target_photons, pulses_target_photon = mt2.NCorrRun("SP", target).photons_singles_tree

theta_bins = [0, 30,33, 54, 67, 77,86, 97,103, 114,127, 135,151, 180]


hist_Cf = TH1F(binarrays=theta_bins)
hist_target = TH1F(binarrays=theta_bins)
hist_Al = TH1F(binarrays=theta_bins)
hist_target_photons = TH1F(binarrays=theta_bins)
hist_target_noise = TH1F(binarrays=theta_bins)
hist_Al_photons = TH1F(binarrays=theta_bins)
hist_target_noise_Al_sub = TH1F(binarrays=theta_bins)

def get_det_photon_rates(tree):
    hits = {angle: 0 for angle in mt.angles}
    hits[-30] = hits[-330] = 0
    pulse_number = 0

    for index in np.linspace(1, tree.GetEntries(), 10000):
        tree.GetEntry(int(index))
        evt = tree
        for hit in evt.photons.hits:

            if hit.ForwardDet != 0:
                det = hit.ForwardDet * hit.det
            else:
                det = hit.det
            hits[det] += 1
        pulse_number = evt.PulseNumber

    for det, n_events in hits.iteritems():
        r = (float(n_events) / pulse_number) * (tree.GetEntries() / 10000.)
        hits[det] = round(r,3)

    return hits

def get_weight_cut(det, w, draw_par):
    draw_par={1:"neutrons", 2:"photons", 3:"noise"}[draw_par]
    if abs(det) in [30, 330]:
        return "{0}*({draw_par}.hits.det == {1} && {draw_par}.hits.ForwardDet == {f})"\
            .format(w, abs(det), draw_par=draw_par, f= -1 if det<0 else 1)
    else:
        return "{0}*({draw_par}.hits.det == {1})".format(w, det, draw_par=draw_par)

def get_DT_cut(tree, draw_par):
    cuts = []
    for det, rate in get_det_photon_rates(tree).iteritems():
        cut = get_weight_cut(det, 1./(1.-rate), draw_par)
        cuts.append(cut)
    return " + ".join(cuts)


Al_subtract_weight = []
noise_AL_subtract_weight = []
for (detAl, rateAl), (det_target, rate_target) in zip(get_det_photon_rates(tree_Al_photons).iteritems(),get_det_photon_rates(tree_target_photons).iteritems()):
    rate_scale = np.log(1./(1-rate_target))/np.log(1./(1-rateAl))
    weight_AL = round(rate_scale*1.0/(1-rateAl),3)

    weight_noise = round(rate_scale * 1.0 / (1 - rate_target),3)

    Al_subtract_weight.append(get_weight_cut(detAl, weight_AL, 1))
    noise_AL_subtract_weight.append(get_weight_cut(detAl, weight_noise, 3))

Al_subtract_weight = " + ".join(Al_subtract_weight)
noise_AL_subtract_weight = " + ".join(noise_AL_subtract_weight)
print(Al_subtract_weight)
print(noise_AL_subtract_weight)


drw = "neutrons.hits.theta_abs*180/3.14"

hist_Cf.Project(treeCf, drw)
hist_target.Project(tree_target, drw, cut="{0} && {1}".format(mt2.get_good_run_cut(target), get_DT_cut(tree_target_photons, 1)), weight=1.0/mt2.get_pulses(mt2.get_good_run_cut(target),target))
hist_Al.Project(tree_Al, drw, cut=Al_subtract_weight ,weight=1.0/pulses_Al)
hist_target_noise_Al_sub.Project(tree_target_noise, "noise.hits.theta_abs*180/3.14", cut="{0} && {1}".format("noise.hits.tof>550",noise_AL_subtract_weight), weight=1.0/pulses_target_noise)
hist_target_noise.Project(tree_target_noise,  "noise.hits.theta_abs*180/3.14", cut="{0} && {1}".format("noise.hits.tof>550",get_DT_cut(tree_target_photons, 3)) ,weight=1.0/pulses_target_noise)

hist_target_noise *= (120.)/(750-550)
hist_target_noise_Al_sub *= (120.)/(750-550)

hist_Al -= hist_target_noise_Al_sub
print("noise: ",sum(hist_target_noise.binvalues)/sum(hist_target.binvalues))
print("photons: ",sum((hist_Al).binvalues)/sum(hist_target.binvalues))
#
hist_target = (hist_target - hist_Al- hist_target_noise)
hist_target /= hist_Cf

hist_target.set_max_rel_err(0.15)

x = np.array([32, 55.7, 78.42, 102, 124, 147], dtype=np.float32)
y = np.array(hist_target.binvalues[np.where(hist_target.binvalues!=0)], dtype=np.float32)
erry = np.array(hist_target.binerrors[np.where(hist_target.binvalues!=0)], dtype=np.float32)

# y[0] *= 0.9

gr = ROOT.TGraphErrors(len(x), x,y,np.zeros_like(x), erry)

gr.SetMarkerStyle(32)
gr.GetXaxis().SetTitle("#theta_{abs}")
gr.GetYaxis().SetTitle("Y(#theta_{{abs}})_{{ {0} }}/Y(#theta_{{abs}})_{{Cf252}}".format(target))
hist_target.SetStats(0)


c1 = ROOT.TCanvas()
mt2.thesis_plot([gr], big_font=0.05)

gr.Draw()
gr.SetTitle("")




if __name__ == "__main__":
    import ROOT as ROOT
    from multiprocessing import Process, Queue
    import time, sys, os
    def input_thread(q, stdin):
        while True:
            print 'ROOT: '
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
                    if ___cmd___ == 'tb':
                        ___tb___ = ROOT.TBrowser()
                    else:
                        try:
                            exec (___cmd___, globals())
                        except:
                            print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()

