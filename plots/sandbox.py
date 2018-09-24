import ROOT
import numpy as np
import mytools as mt

import mytools2 as mt2

from TH1Wrapper import TH1F

target = "DU"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False,
                                                 Forward=True).neutrons_doubles_tree
treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False,
                                                 Forward=True).neutrons_doubles_tree

binwidth = 15

cuts = ['', 'neutrons.coinc_hits.ForwardDet ==0']
histos = []

for i, cut in enumerate(cuts):
    histSP = TH1F(20, 180, binwidths=binwidth)
    histDP = TH1F(20, 180, binwidths=binwidth)

    histos.append(histSP)

    evnts = histSP.Project(treeSP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', cut,
                           weight=1.0 / pulses_SP_doubles)
    histDP.Project(treeDP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', cut, weight=1.0 / pulses_DP_doubles)

    print(cut, evnts)

    histSP -= 0.5 * histDP

    histSP /= (0.5 * histDP)

    histSP.Draw('' if i == 0 else 'same')

    histSP.SetMinimum(0)

    if i == 0:
        histSP.SetLineColor(ROOT.kRed)


for i, (x1,x2) in enumerate(zip(histos[0].binvalues,histos[1].binvalues)):
    print('Diff {0}: {1:.2f} sig = {2:.2f}'.format(histSP.bincenters[0][i],100*(x1-x2)/x1, (x1-x2)/histSP.binerrors[i]))

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
#
