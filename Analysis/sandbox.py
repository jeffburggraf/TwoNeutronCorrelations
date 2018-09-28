import ROOT
from TH1Wrapper import TH1F, TH2F
import numpy as np
import mytools2 as mt2
import mytools as mt



file = ROOT.TFile("test.root","recreate")

a = np.array([0.0], np.float32)
tree = ROOT.TTree("tree", "tree")
tree.Branch("a", a, "a/F")


for i in np.random.uniform(0,100,10000):
    a[0] = i
    tree.Fill()

hist = TH1F(0,100,binwidths=1)



l = np.linspace(0,100,10)
cuts = []

for i1,i2 in zip(l[:-1],l[1:]):
    print(i1,i2)
    cuts.append("{0}*({1})".format((i1+i2)/2.,mt.cut_rangeAND([i1,i2], "a")))

cut = " + ".join(cuts)

print(cut)

hist.Project(tree, "a", cut = cut, weight=1.0/10000)
tb = ROOT.TBrowser()

hist.Draw()





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




