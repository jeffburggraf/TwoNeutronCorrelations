import ROOT
import numpy as np
import mytools2 as mt2
from TH1Wrapper import TH1F
import mytools as mt
from matplotlib import pyplot as plt

treeSP, _ = mt2.NCorrRun("SP","Cf252",generate_dictionary=False,Forward = False).all_singles_tree

# fig,ax = plt.subplots(4,3)


for i, det in enumerate(mt.angles):
    hist = TH1F(-30, 40, binwidths=0.3)

    cut = "trigger_array[0]!=0 && trigger_array[1]!=0 && hits.phi == {0}".format(det)
    hist.Project(treeSP, "0.5*(all.hits.top + all.hits.bot - trigger_array[0] - trigger_array[1])",cut)

    # plt.subplot(4,3,i+1)

    plt.plot(hist.bincenters[0], hist.binvalues, label=det)
    print('wtf')

plt.legend()
# # hist.Project(treSP, 'trig[1] - trig[2]')
# # hist.Draw()
#
# print("dfij")
# # tb = ROOT.TBrowser()
#



plt.ion()
plt.show()

# print("dkdk")


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