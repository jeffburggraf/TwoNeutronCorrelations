import ROOT
import numpy as np
import os
import mytools as mt
# =================
min_erg = 0.4
max_events = None
# =================

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


neutron_struct = ROOT.neutrons()
PulseNumber = np.array([0], dtype=np.int32)

if not os.path.isdir("/Volumes/JeffMacEx/2nCorrData/Production/Forward/SP_FREYA"):
    os.mkdir("/Volumes/JeffMacEx/2nCorrData/Production/Forward/SP_FREYA/")
if not os.path.isdir("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA"):
    os.mkdir("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/")

_F = ROOT.TFile("/Users/jeffreyburggraf/FREYA/MyFREYA/FREYA.root")
tree_FREYA = _F.Get("tree")

fSP = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/SP_FREYA/SP_FREYA_neutrons_coinc.root", 'RECREATE')
treeSP = ROOT.TTree('tree', "")
treeSP.Branch('neutrons', neutron_struct)
treeSP.Branch('PulseNumber', PulseNumber, "PulseNumber/I")

def runSP():
    eventsSP = 0
    fSP.cd()
    fissID = 0

    for FREYA_evt in tree_FREYA:
        PulseNumber[0] = fissID
        fissID += 1
        neutron_struct.xs = FREYA_evt.inducing_par.xs

        coinc_hit = ROOT.coinc_hit()
        nu = FREYA_evt.neutronNu

        if nu<2:continue
        _a = np.arange(nu)
        np.random.shuffle(_a)
        _a = _a[:len(_a) - len(_a)%2]

        indicies = np.split(_a, nu//2)
        for Is in indicies:
            i1, i2 = tuple(Is)
            coinc_hit.erg[0] = FREYA_evt.neutronEnergies[i1]
            coinc_hit.erg[1] = FREYA_evt.neutronEnergies[i2]
            if coinc_hit.erg[0]< min_erg or coinc_hit.erg[1]< min_erg:continue
            xyz1 = [FREYA_evt.neutronDirx[i1], FREYA_evt.neutronDiry[i1], FREYA_evt.neutronDirz[i1]]
            xyz2 = [FREYA_evt.neutronDirx[i2], FREYA_evt.neutronDiry[i2], FREYA_evt.neutronDirz[i2]]
            coinc_hit.coinc_theta = mt.vectorAngle(xyz1, xyz2)
            coinc_hit.theta_abs[0] = mt.vectorAngle([1,0,0], xyz1)
            coinc_hit.theta_abs[1] = mt.vectorAngle([1,0,0], xyz2)
            neutron_struct.coinc_hits.push_back(coinc_hit)
            eventsSP +=1

            treeSP.Fill()

            if max_events and eventsSP>max_events:
                treeSP.Write()
                neutron_struct.coinc_hits.clear()

                return


        neutron_struct.coinc_hits.clear()


    treeSP.Write()
runSP()
fSP.Close()

fDP = ROOT.TFile("/Volumes/JeffMacEx/2nCorrData/Production/Forward/DP_FREYA/DP_FREYA_neutrons_coinc.root", 'RECREATE')
treeDP = ROOT.TTree('tree', "")
treeDP.Branch('neutrons', neutron_struct)
treeDP.Branch('PulseNumber', PulseNumber, "PulseNumber/I")


def runDP():
    fDP.cd()
    eventsDP = 0
    entry = 1
    fissID = 0
    while entry < tree_FREYA.GetEntries():
        PulseNumber[0] = fissID
        fissID += 1
        i1, i2 = entry, entry+1
        tree_FREYA.GetEntry(i1)
        xs1 = tree_FREYA.inducing_par.xs
        nu1 = tree_FREYA.neutronNu
        shuffled_indicies1 = np.arange(nu1)
        np.random.shuffle(shuffled_indicies1)
        xyzs1 = [[tree_FREYA.neutronDirx[i],tree_FREYA.neutronDiry[i],tree_FREYA.neutronDirz[i]] for i in shuffled_indicies1]
        ergs1 = [tree_FREYA.neutronEnergies[i] for i in shuffled_indicies1]

        tree_FREYA.GetEntry(i2)
        xs2 = tree_FREYA.inducing_par.xs
        nu2 = tree_FREYA.neutronNu
        shuffled_indicies2 = np.arange(nu2)
        np.random.shuffle(shuffled_indicies2)
        xyzs2 = [[tree_FREYA.neutronDirx[i],tree_FREYA.neutronDiry[i], tree_FREYA.neutronDirz[i]] for i in shuffled_indicies2]
        ergs2 = [tree_FREYA.neutronEnergies[i] for i in shuffled_indicies2]

        neutron_struct.xs = xs1*xs2
        entry += 2

        if nu2 < nu1:
            xyzs2, xyzs1 = xyzs1, xyzs2
            ergs2, ergs1 = ergs1 ,ergs2
            nu2,     nu1 =   nu1,   nu2

        if not (nu1 > 0 and nu2 > 0):
            continue

        for index in range(nu1):
            coinc_hit = ROOT.coinc_hit()
            coinc_hit.erg[0] = ergs1[index]
            coinc_hit.erg[1] = ergs2[index]
            if coinc_hit.erg[0]< min_erg or coinc_hit.erg[1]< min_erg:continue

            coinc_hit.theta_abs[0] = mt.vectorAngle(xyzs1[index], [1,0,0])
            coinc_hit.theta_abs[1] = mt.vectorAngle(xyzs2[index], [1,0,0])
            coinc_hit.coinc_theta = mt.vectorAngle(xyzs1[index], xyzs2[index])
            neutron_struct.coinc_hits.push_back(coinc_hit)
            eventsDP += 1

        treeDP.Fill()
        neutron_struct.coinc_hits.clear()
        if max_events and eventsDP>max_events:
            treeDP.Write()
            return

    treeDP.Write()

runDP()
fDP.Close()

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