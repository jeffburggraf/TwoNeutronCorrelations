import mytools2 as mt2
from TH1Wrapper import TH1F
import mytools as mt
import ROOT

n = 10000

tree, _ = mt2.NCorrRun('SP', "DU").photons_singles_tree


tree.GetEntry(n)
n_pulses = tree.PulseNumber

histos = []


tot = 0
tot2 = 0
for angle in mt.angles:
    cut = "photons.hits.det == {}".format(angle)
    if abs(angle) in [30, 330]:
        if angle<0:
            cut += " && photons.hits.ForwardDet == -1"
        else:
            cut += " && photons.hits.ForwardDet == 1"

    events = tree.Project("", "", cut, "", n)

    if abs(angle) in [30, 330]:
        e = events / float(n) * 0.333
    else:
        e = events / float(n)


    tot += (1-e)
    tot2 += e

print tot

for i in range(10000):
    tree.GetyEmntry()
    if tree.k:0
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

