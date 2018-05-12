import ROOT
import numpy as np
import mytools2 as mt2

target = "D2O"

treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

treeSP_singles, pulses_SP_singles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_singles_tree

treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree


zero_hits = (pulses_SP_doubles-treeSP_singles.GetEntries())/pulses_SP_doubles
zero_hits_err = np.sqrt(pulses_SP_doubles-treeSP_singles.GetEntries())/pulses_SP_doubles

single_hits = (treeSP_singles.GetEntries() - treeSP_doubles.GetEntries())/pulses_SP_doubles
single_hits_err = np.sqrt(treeSP_singles.GetEntries() - treeSP_doubles.GetEntries())/pulses_SP_doubles

doubles_hits = treeSP_doubles.GetEntries() / pulses_SP_doubles
doubles_hits_err = np.sqrt(treeSP_doubles.GetEntries()) / pulses_SP_doubles

triple_hits = treeSP_doubles.Project("","","neutrons.nhits ==3 ") / pulses_SP_doubles
triple_hits_err = np.sqrt(treeSP_doubles.Project("","", "neutrons.nhits ==3 ")) / pulses_SP_doubles

x = np.array([0,1,2],np.float64)
x_err = np.zeros_like(x)

y = np.array([zero_hits, single_hits, doubles_hits], np.float64)
y_err = 1.5*np.array([zero_hits_err, single_hits_err, doubles_hits_err],np.float64)

func = ROOT.TF1("fit", "TMath::Poisson(x, [0])")
# func = ROOT.TF1("fit", "[0]*2**([1]*x)/x")
func.SetParameter(0,single_hits)
func.SetParameter(1,0)
func.SetParName(0,"#lambda")

gr = ROOT.TGraphErrors(len(x), x, y, x_err, y_err)
gr.Draw("*A")
gr.Fit("fit")

gr.GetYaxis().SetRangeUser(y[-1]*0.5,1.1);
gr.GetXaxis().SetLimits(-0.2,3);

gr.GetXaxis().SetTitle("Number of neutrons detected in coincidence");
gr.GetYaxis().SetTitle("Counts per pulse");

leg = ROOT.TLegend()
leg.AddEntry(func, "Poissonian fit")
leg.AddEntry(gr, "D2O data")

leg.Draw()

gr.SetMarkerStyle(32)
ROOT.gPad.Update()

mt2.thesis_plot(gr)

ROOT.gPad.SetLogy();

ROOT.gStyle.SetOptFit(1111)

# TB = ROOT.TBrowser()


print y,y_err


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