import ROOT
import numpy as np
import mytools2 as mt2
import itertools
import mytools as mt

# =-------------------
target = "DU"
max_events = None
sliding_window_width = 40
# =-------------------
np.random.seed(1)

ROOT.gROOT.ProcessLine(".L twontree.h")

old_tree, _=mt2.NCorrRun("SP",target).neutrons_doubles_tree


class Pulse:
    def __init__(self):
        self.ergs = [old_tree.neutrons.coinc_hits[0].erg[0], old_tree.neutrons.coinc_hits[0].erg[1]]
        self.xs = [old_tree.neutrons.coinc_hits[0].x[0], old_tree.neutrons.coinc_hits[0].x[1]]
        self.ys = [old_tree.neutrons.coinc_hits[0].y[0], old_tree.neutrons.coinc_hits[0].y[1]]
        self.zs = [old_tree.neutrons.coinc_hits[0].z[0], old_tree.neutrons.coinc_hits[0].z[1]]
        self.theta_abs = [old_tree.neutrons.coinc_hits[0].theta_abs[0], old_tree.neutrons.coinc_hits[0].theta_abs[1]]
        self.det = [old_tree.neutrons.coinc_hits[0].det[0], old_tree.neutrons.coinc_hits[0].det[1]]
        self.PulseNumber =old_tree.PulseNumber

def get_uncorr_events(pulse1, pulse2):
    assert isinstance(pulse2, Pulse)
    assert isinstance(pulse1, Pulse)
    result = {"n_events":0,"theta_nn":[],"ergs":[], "theta_abs":[]}
    for i in range(2):
        for j in range(2):
            if pulse1.det[i] == pulse2.det[j]:
                continue

            result["n_events"] += 1
            result["ergs"].append([0,0])
            result["theta_abs"].append([0,0])

            x1,y1,z1 = pulse1.xs[i],pulse1.ys[i],pulse1.zs[i]
            x2,y2,z2 = pulse2.xs[j],pulse2.ys[j],pulse2.zs[j]

            result["theta_nn"].append(mt.vectorAngle((x1,y1,z1),(x2,y2,z2), radians=False))
            if result["theta_nn"][-1]==0:
                print (x1,y1,z1), (x2,y2,z2), pulse1.det, pulse2.det

            result["ergs"][-1][0] = pulse1.ergs[0]
            result["ergs"][-1][1] = pulse2.ergs[1]

            result["theta_abs"][-1][0] = pulse1.theta_abs[0]
            result["theta_abs"][-1][1] = pulse2.theta_abs[1]

    result["diff"] = abs(pulse1.PulseNumber - pulse2.PulseNumber)/240.

    return result

pules = []
for i in range(old_tree.GetEntries()):
    old_tree.GetEntry(i+1)
    pules.append(Pulse())

w = int(sliding_window_width/2.)

new_file = ROOT.TFile(
    "/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/New_trees/{}.root".format(target),
    "recreate")
new_tree = ROOT.TTree("tree", "tree")

br_erg = np.zeros(2, dtype=np.float32)
br_theta_abs = np.zeros(2, dtype=np.float32)
br_theta_nn = np.zeros(1, dtype=np.float32)
br_diff = np.zeros(1, dtype=np.float32)

new_tree.Branch("theta_nn",br_theta_nn, "theta_nn/F")
new_tree.Branch("abs_theta",br_theta_abs, "abs_theta[2]/F")
new_tree.Branch("erg",br_erg, "erg[2]/F")
new_tree.Branch("time_diff",br_diff, "time_diff/F")


print("entries: {0}".format(old_tree.GetEntries()))


def loop():
    n_events = 0
    next_print = 10000

    for i_center in range(w, old_tree.GetEntries()-w,w):
        for i1 in range(i_center-w, i_center+w-1):
            for i2 in range(i1, i_center+w):
                pulse1 = pules[i1]
                pulse2 = pules[i2]
                result = get_uncorr_events(pulse1, pulse2)
                br_diff[0]=result["diff"]

                for i in xrange(result["n_events"]):
                    br_theta_nn[0] = result['theta_nn'][i]
                    br_erg[0] = result['ergs'][i][0]
                    br_erg[1] = result['ergs'][i][1]

                    br_theta_abs[0] = result['theta_abs'][i][0]
                    br_theta_abs[1] = result['theta_abs'][i][1]
                    new_tree.Fill()

                n_events += 4

                if n_events>next_print:
                    print("{:.1f}% complete...".format(100*float(i_center)/old_tree.GetEntries()))
                    next_print +=10000

                if max_events is not None and n_events>max_events:return
loop()
new_tree.Write()

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

