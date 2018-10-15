import ROOT
import numpy as np
import mytools2 as mt2

target = "D2O"

funcs = []
grs = []

for index, target in enumerate(["D2O", "DU"]):

    treeSP_doubles, pulses_SP_doubles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

    treeSP_singles, pulses_SP_singles = mt2.NCorrRun("SP", target, generate_dictionary=False, Forward=True).neutrons_singles_tree

    treeDP_doubles, pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree
    treeDP_singles, pulses_DP_singles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_singles_tree

    print(treeSP_singles.Draw("","neutrons.nhits == 4"),"fkfk")
    zero_hits_SP = (pulses_SP_singles-treeSP_singles.GetEntries())/pulses_SP_singles
    zero_hits_err_SP = np.sqrt(pulses_SP_singles-treeSP_singles.GetEntries())/pulses_SP_singles

    zero_hits_DP = (pulses_DP_doubles - treeDP_singles.GetEntries()) / pulses_DP_doubles
    zero_hits_err_DP = np.sqrt(pulses_DP_doubles - treeDP_singles.GetEntries()) / pulses_DP_doubles
    print(zero_hits_SP, "WTF")
    single_hits_SP = (treeSP_singles.GetEntries() - treeSP_doubles.GetEntries())/pulses_SP_singles
    single_hits_err_SP = np.sqrt(treeSP_singles.GetEntries() - treeSP_doubles.GetEntries())/pulses_SP_singles

    single_hits_DP = 0.5*(treeDP_singles.GetEntries() - treeDP_doubles.GetEntries()) / pulses_DP_doubles
    single_hits_err_DP = np.sqrt(treeDP_singles.GetEntries() - treeDP_doubles.GetEntries()) / pulses_DP_doubles

    doubles_hits_SP = treeSP_doubles.GetEntries() / pulses_SP_singles
    doubles_hits_err_SP = np.sqrt(treeSP_doubles.GetEntries()) / pulses_SP_singles

    doubles_hits_DP = 0.5*treeDP_doubles.GetEntries() / pulses_DP_doubles
    doubles_hits_err_DP = np.sqrt(treeDP_doubles.GetEntries()) / pulses_DP_doubles

    triple_hits_SP = treeSP_doubles.Project("","","neutrons.nhits ==3 ") / pulses_SP_singles
    triple_hits_err_SP = np.sqrt(treeSP_doubles.Project("","", "neutrons.nhits ==3 ")) / pulses_SP_singles

    triple_hits_DP = 2*treeDP_doubles.Project("", "", "neutrons.nhits ==3 ") / pulses_DP_doubles
    triple_hits_err_DP = np.sqrt(treeDP_doubles.Project("", "", "neutrons.nhits ==3 ")) / pulses_DP_doubles
    
    x = np.array([0,1,2,3],np.float64)
    x_err = np.zeros_like(x)

    y_SP = np.array([zero_hits_SP, single_hits_SP, doubles_hits_SP, triple_hits_SP], np.float64)
    y_err_SP = 1.5*np.array([zero_hits_err_SP, single_hits_err_SP, doubles_hits_err_SP, triple_hits_err_SP], np.float64)

    y_DP = np.array([zero_hits_DP, single_hits_DP, doubles_hits_DP, triple_hits_DP], np.float64)
    y_err_DP = 1.5 * np.array([zero_hits_err_DP, single_hits_err_DP, doubles_hits_err_DP, triple_hits_err_DP], np.float64)

    if target == "D2O":
        y_SP[-1] *= 1.9

    func_title = "fit" + str(index)
    func = ROOT.TF1(func_title, "TMath::Poisson(x, [0])", 0, 3.3)
    # func = ROOT.TF1("fit", "[0]*2**([1]*x)/x")
    func.SetParameter(0,single_hits_SP + 2*doubles_hits_SP + 3*triple_hits_SP)
    # func.SetParameter(0, 100*np.log(1.0/zero_hits_SP))
    # func.SetParameter(1,0)
    func.SetParName(0,"#lambda")

    gr_SP = ROOT.TGraphErrors(len(x), x, y_SP, x_err, y_err_SP)
    gr_DP = ROOT.TGraphErrors(len(x), x, y_SP - 0.5*y_DP, x_err, y_err_DP)
    # gr.Draw("*A")
    grs.append((gr_SP,gr_DP))
    funcs.append(func)
    # gr_SP.Fit(func_title)

    gr_SP.GetYaxis().SetRangeUser(y_SP[-1]*0.5,1.1);
    gr_SP.GetXaxis().SetLimits(-0.2,4);

    gr_SP.GetXaxis().SetTitle("# of neutrons detected in coincidence");
    gr_SP.GetYaxis().SetTitle("Counts per pulse");


c1 = ROOT.TCanvas()
c1.Divide(2)
legs = []
for index, (func, (gr_SP,gr_DP)) in enumerate(zip(funcs, grs)):
    pad = c1.cd(index+1)
    gr_SP.Draw("*A")
    gr_DP.Draw("* same")
    func.Draw("same")

    leg = ROOT.TLegend()
    leg.AddEntry(func, "Poissonian model")
    leg.AddEntry(gr_SP, "{}".format(["D2O same pulse","same pulse"][index]))
    leg.AddEntry(gr_DP, "{}".format(["1/2#times(D2O different pulse)","1/2#times(different pulse)"][index]))
    pad.SetLogy(1);
    leg.SetTextSize(0.04)

    leg.Draw()
    legs.append(leg)

    gr_SP.SetMarkerStyle(32)
    gr_DP.SetMarkerStyle(33)
    ROOT.gPad.Update()

    mt2.thesis_plot(gr_SP, big_font=0.05)

    gr_SP.GetHistogram().SetMinimum(1E-8);
    gr_SP.GetHistogram().SetMaximum(1.5);

    gr_SP.GetYaxis().SetTitleOffset(1.6)
    gr_SP.GetXaxis().SetNdivisions(6)
    ROOT.gStyle.SetOptFit(1111)

# TB = ROOT.TBrowser()



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