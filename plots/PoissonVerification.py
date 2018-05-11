import ROOT
import numpy as np
import mytools2 as mt2



treeSP_doubles, _ = mt2.NCorrRun("SP", "DU", generate_dictionary=False, Forward=True).neutrons_doubles_tree

treeSP_singles, pulses_SP_singles = mt2.NCorrRun("SP", "DU", generate_dictionary=False, Forward=True).neutrons_singles_tree

treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", "DU", generate_dictionary=False, Forward=True).neutrons_doubles_tree


zero_hits = (pulses_SP_singles-treeSP_singles.GetEntries() - treeSP_doubles.GetEntries())/pulses_SP_singles
zero_hits_err = np.sqrt(pulses_SP_singles-treeSP_singles.GetEntries() - treeSP_doubles.GetEntries())/pulses_SP_singles

single_hits = (treeSP_singles.GetEntries())/pulses_SP_singles
single_hits_err = np.sqrt(treeSP_singles.GetEntries())/pulses_SP_singles

doubles_hits = treeSP_doubles.GetEntries() / pulses_DP_doubles
doubles_hits_err = np.sqrt(treeSP_doubles.GetEntries()) / pulses_DP_doubles

print (zero_hits, single_hits, doubles_hits)


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