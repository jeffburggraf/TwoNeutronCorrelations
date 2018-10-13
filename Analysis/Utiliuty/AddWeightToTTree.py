import ROOT
import os
from TH1Wrapper import TH2F, TH1F
import re
import mytools2 as mt2
import time
import numpy as np

# dir = "/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/Analysis/Other"
dir = "/Volumes/JeffMacEx/2nCorrData/Production/Forward/"

ROOT.gROOT.ProcessLine(".L twontree.h+")

ergs = TH1F(0.4, 12, nbinss=20)
histos = []

def gen_FREYA_dict():
    ROOT.gROOT.ProcessLine(
    """
    struct coinc_hit {
        float erg[2]; 
        double theta_abs[2];
        double coinc_theta;
    };
    struct neutrons {
        std::vector<coinc_hit> coinc_hits;
        double xs;
    };
    
    struct FF{
        int A,Z;
        double KE, dir[3];
    };
    
    struct inducing_par{
          double erg;
          double dir[3];
          double xs;
          int type;
    };
    
    """)


def closure(target):
    freya = True if "FREYA" in target else False
    if freya:
        gen_FREYA_dict()
    subtract_accidentals = False if (freya or target == "Cf252") else True

    __file_SP = ROOT.TFile(__file_path_SP)
    __file_DP = ROOT.TFile(__file_path_DP)

    treeSP = __file_SP.Get("tree")
    treeDP = __file_DP.Get("tree")

    treeSP.GetEntry(treeSP.GetEntries() - 1)
    n_pulsesSP = float(treeSP.PulseNumber)

    treeDP.GetEntry(treeDP.GetEntries() - 1)
    n_pulsesDP = float(treeDP.PulseNumber)

    n_bins = 12

    # while n_bins>0:
    # __erg_hist__SP = TH1F(0.4, 9, 100)
    # __erg_hist__SP.Project(treeSP, "neutrons.coinc_hits.erg",  weight=1.0/n_pulsesSP)
    #
    # if subtract_accidentals:
    #     __erg_hist__DP = TH1F(0.4, 9, 100)
    #     __erg_hist__DP.Project(treeDP, "neutrons.coinc_hits[].erg", weight=1.0/n_pulsesDP)
    #     __erg_hist__SP -= 0.5*__erg_hist__DP
    #
    # erg_bins, _ = mt2.median(__erg_hist__SP, n_bins)
    erg_bins = np.concatenate((np.linspace(0, 5, n_bins), np.linspace(5, 9, 4)))

    histSP = TH2F(binarrays=erg_bins)
    histDP = TH2F(binarrays=erg_bins)

    histSP.Project(treeSP, "neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]", weight=1.0 / n_pulsesSP)

    histDP.Project(treeDP, "neutrons.coinc_hits[0].erg[0]:neutrons.coinc_hits[0].erg[1]", weight=1.0 / n_pulsesDP)

    if subtract_accidentals:
        histSP -= histDP * 0.5

    histSP = 0.5 * (histSP + histSP.transpose())
    histDP = 0.5 * (histDP + histDP.transpose())

    weighted_hist = histSP / (0.5 * histDP)

    #     MAX_REL_ERROR = max((weighted_hist.binerrors / weighted_hist.binvalues).flatten())
    #
    #     if MAX_REL_ERROR<0.35:
    #         break
    #     else:
    #         n_bins -= 1
    #
    # if n_bins != 7:
    #     print("Lowering the number of bins for {target} from 7 to {n_bins} because max_rel_error was {0:.1f}". \
    #         format(100 * MAX_REL_ERROR, target=target, n_bins=n_bins))
    # if n_bins == 1:
    #     weighted_hist = TH2F(binarrays=erg_bins)
    #     weighted_hist += 1
    # else:
    #     print("\nMax rel error for {target} is  {0:.1f}%".format(100*MAX_REL_ERROR, target=target))

    weighted_hist.SetTitle(target)

    weighted_hist.Smooth()
    weighted_hist.update_bin_containers_from_hist()

    # weighted_hist.DrawErrrorBars()

    weighted_hist.SetMinimum(0)

    __new_file__ = ROOT.TFile("{0}/{1}/DP_weights.root".format(dir, dir_name), "recreate")
    weight_tree = ROOT.TTree("tree", "")

    weight_br = np.array([0], dtype=np.float32)
    weight_tree.Branch("DP_weight", weight_br, "DP_weight/F")

    n_events = 0
    max_events = None

    evt = treeDP
    zer_lengths = 0
    overflows = 0

    original_weight = 0

    weights = np.zeros(treeDP.GetEntries(), dtype=np.float32)
    for tree_i in range(treeDP.GetEntries()):
        treeDP.GetEntry(tree_i)
        if len(evt.neutrons.coinc_hits) == 0:
            weights[tree_i] = 0
            zer_lengths += 1
            continue
        erg_1 = evt.neutrons.coinc_hits[0].erg[0]
        erg_2 = evt.neutrons.coinc_hits[0].erg[1]

        n_events += 1

        i1 = i2 = None
        for i, b in enumerate(weighted_hist.__binRightEdges__[0]):
            if i1 is not None and i2 is not None: break
            if i1 is None and b > erg_1:
                i1 = i
            if i2 is None and b > erg_2:
                i2 = i
        if i1 is None:
            overflows += 1
        if i2 is None:
            overflows += 1

        if (i1 is None) or (i2 is None):
            weights[tree_i] = 0
            continue

        weights[tree_i] = weighted_hist.binvalues[i1][i2]
        original_weight += 1

        if max_events is not None and n_events > max_events:
            break

    weights *= original_weight / np.sum(weights)
    for w in weights:
        weight_br[0] = w
        weight_tree.Fill()

    weight_tree.Write()

    histos.append(weighted_hist)
    weighted_hist.Draw("lego hist E")
    weighted_hist.DrawErrrorBars()
    weighted_hist.SetTitle(target)

    return weighted_hist


___g___ = ROOT.gSystem.ProcessEvents
for dir_name in os.listdir(dir):
    _m = re.match("DP_(.+)", dir_name)

    if _m:
        target = _m.groups(1)[0]
        __file_path_DP = "{dir}/{dir_name}/DP_{target}_neutrons_coinc.root".format(dir_name=dir_name, dir=dir, target=target)
        if os.path.exists("{dir}/{dir_name}/DP_weights.root".format(dir_name=dir_name, dir=dir, target=target)):
            print("'{dir_name}/DP_weights.root' already exists!".format(dir_name=dir_name, dir=dir, target=target))
            continue

        __file_path_SP = "{dir}/SP_{target}/SP_{target}_neutrons_coinc.root".format(dir_name=dir_name, dir=dir, target=target)
        if not os.path.exists(__file_path_DP) or not os.path.exists(__file_path_SP):
            print ("Coincident neutron trees for {} does not exist!".format(target))
            continue
    else:
        continue
    print("continuing for target {}".format(target))
    try:
        hist = closure(target)
    except Exception as e:
        os.remove("{dir}/{dir_name}/DP_weights.root".format(dir_name=dir_name, dir=dir, target=target))
        print(e)
        continue


    # hist.Draw("lego hist E")
    # ___g___()



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
