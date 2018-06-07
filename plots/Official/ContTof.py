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


dx = (2.54*15)/3

ps =  [-3*dx, -2*dx,-dx,0,dx,2*dx, 3*dx]


x_mean = []
y_mean = []
erry_mean = []

y_top = []
y_bot = []
erry_top = []
erry_bot = []

offset = 2.54*4
for p in ps:
    d1 = abs(15*2.54 - p)
    d2 = abs(-15*2.54 - p)
    hist = TH1F(-10, 10, binwidths=0.25)

    mean_times =[]
    tops = []
    bots = []
    for delta in range(70):
        t1 = d1/(1.52E10)*1E9 + np.random.randn()*1.5
        t2 = d2/(1.52E10)*1E9 + np.random.randn()*1.5

        Dt = t2+t1
        mean_times.append(0.5*Dt)
        tops.append(t1)
        bots.append(t2)

    x_mean.append(p)#- (2.54*30 + offset)/2.)

    y_mean.append(np.mean(mean_times))
    y_top.append(np.mean(tops))
    y_bot.append(np.mean(bots))

    erry_mean.append(np.std(mean_times)/4.)
    erry_top.append(np.std(tops)/4.)
    erry_bot.append(np.std(bots)/4.)

x_mean = np.array(x_mean, dtype=np.float64)
y_mean = np.array(y_mean, dtype=np.float64)
erry_mean = np.array(erry_mean, dtype=np.float64)

y_top = np.array(y_top, dtype=np.float64)
y_bot = np.array(y_bot, dtype=np.float64)
erry_top = np.array(erry_top, dtype=np.float64)
erry_bot = np.array(erry_bot, dtype=np.float64)

ot = 10,15

y_mean += np.mean(ot)
y_bot += ot[0]
y_top+= ot[1]

gr = ROOT.TGraphErrors(len(x_mean), x_mean, y_mean,np.ones_like(erry_mean)*0, erry_mean)
gr_top = ROOT.TGraphErrors(len(x_mean), x_mean, y_top,np.ones_like(erry_mean)*0, erry_top)
gr_bot = ROOT.TGraphErrors(len(x_mean), x_mean, y_bot, np.ones_like(erry_mean)*0, erry_bot)

gr.Draw('A*')
gr_bot.Draw('*same')
gr_top.Draw('*same')

gr.SetMinimum(0)
# gr.SetRange(0)

gr.GetXaxis().SetTitle('Distance from detector center [cm]')
gr.GetYaxis().SetTitle('PMT timing [ns]')
gr.SetMarkerStyle(33)
gr_top.SetMarkerStyle(26)
gr_bot.SetMarkerStyle(32)


mt2.thesis_plot(gr, True)

# leg = ROOT.TLegend()
# leg.AddEntry(gr, ' ','ep')
# leg.AddEntry(gr_top, ' ','ep')
# leg.AddEntry(gr_bot, ' ','ep')
# leg.SetTextSize(0.06)
# leg.Draw()
# leg.SetBorderSize(0);

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