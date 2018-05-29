import ROOT
import numpy as np
import mytools2 as mt2
from TH1Wrapper import TH1F
import mytools as mt
import pickle
from sklearn.linear_model import Ridge


targets = "DU", 'Cf252'

treeSP1, pulses_SP1_doubles = mt2.NCorrRun("SP", targets[0], generate_dictionary=False, Forward=True).neutrons_singles_tree
treeSP2, pulses_SP2_doubles = mt2.NCorrRun("SP", targets[-1], generate_dictionary=False,Forward=True).neutrons_singles_tree

with open('Al.pickle', 'r') as file:
    features = pickle.load(file)
print(features)

binsphi = [0]
for degree in mt.angles:
    x1 = degree - 12
    x2 = degree + 12

    if x2>180:
        x2 = 180

    binsphi.append(x1)
    binsphi.append(x2)

    if x2>=180:break

binsphi = set([0] + binsphi)
binsphi = sorted(list(binsphi))
print(binsphi)

erg_hist1 = TH1F(0.35,6,binwidths = 0.5)
erg_hist2 = TH1F(0.35,6,binwidths=0.5)
erg_hist1.Project(treeSP1, 'neutrons.hits.erg')
erg_hist2.Project(treeSP1, 'neutrons.hits.erg')

erg_hist1 += erg_hist2
erg_bins,_ = mt2.median(erg_hist1,3)


c1 = ROOT.TCanvas()
c1.Divide(2,2)

i=1
for E1,E2 in zip(erg_bins[:-1],erg_bins[1:]):
    c1.cd(i)
    hist1 = TH1F(binarrays=binsphi)
    hist2 = TH1F(binarrays=binsphi)
    cut = mt.cut_rangeAND([E1, E2], 'neutrons.hits.erg')
    hist1.Project(treeSP1, '180/3.1415*neutrons.hits.theta_abs', cut)

    hist2.Project(treeSP2, '180/3.1415*neutrons.hits.theta_abs', cut)
    hist1 /= hist2

    hist1.SetTitle('{0}/{1}  {2}<E<{3}'.format(targets[0], targets[1],round(E1, 2),round(E2, 2)))
    hist1.GetXaxis().SetTitle("#theta_{abs}")
    hist1.GetYaxis().SetTitle("ratio")

    hist1.Draw(make_new_canvas=False)
    mt2.thesis_plot(hist1)
    i+=1

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
                    try:
                        exec (___cmd___, globals())
                    except:
                        print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()


