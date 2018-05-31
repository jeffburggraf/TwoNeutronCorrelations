import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2

import ROOT
import numpy as np
import mytools2 as mt2

# target = "cf252"
# treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
# hist = TH1F(-30,30, binwidths=1)
# hist.Project(treeSP_doubles, 'neutrons.pmt_diff', 'neutrons.pmt_diff!=0')
# hist.Draw()
# TB = ROOT.TBrowser()


dx = (2.54*30)/5

x = []
y = []
erry = []
offset = 2.54*4
for p in np.arange(offset,offset + dx*5,dx):
    d1 = 30*2.54 - p
    d2 = p
    hist = TH1F(-10, 10, binwidths=0.25)
    print(p)

    entries =[]
    for delta in range(150):
        t1 = d1/(1.52E10)*1E9 + np.random.randn()*1
        t2 = d2/(1.52E10)*1E9 + np.random.randn()*1

        Dt = t2-t1
        entries.append(Dt)
    x.append(d2 - 40.64)#- (2.54*30 + offset)/2.)
    y.append(np.mean(entries) + 10)
    erry.append(np.std(entries)/np.sqrt(len(entries)))

x = np.array(x, dtype=np.float64)
y = np.array(y, dtype=np.float64)
erry = np.array(erry, dtype=np.float64)

gr = ROOT.TGraphErrors(len(x), x, y,np.ones_like(erry)*0.5, erry)

gr.Draw('A*')
f = gr.Fit('pol1', 'S')

gr.GetXaxis().SetTitle('Distance from detector center [cm]')
gr.GetYaxis().SetTitle('PMT timing difference [ns]')
gr.SetMarkerStyle(31)

# hist.Draw()
mt2.thesis_plot(gr)
leg = ROOT.TLegend()
leg.AddEntry(gr, 'Measurement', 'lep')
leg.AddEntry('d', 'Fit: y = ax + b')
leg.Draw()


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