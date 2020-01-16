import ROOT
import numpy as np
import mytools as mt
import mytools2 as mt2
from TH1Wrapper import TH1F
import matplotlib as mpl
from matplotlib import pyplot as plt
font = {'family': 'DejaVu Sans',
        'size': 16}
## for Palatino and other serif fonts use:

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc("savefig", dpi=300)

target = "Al"

Al_treeDP_doubles, Al_pulses_DP_doubles = mt2.NCorrRun("DP", target, generate_dictionary=False, Forward=True).neutrons_doubles_tree

DU_treeDP_doubles, DU_pulses_DP_doubles = mt2.NCorrRun("DP", 'DU', generate_dictionary=False, Forward=True).neutrons_doubles_tree

erg_bins = [0, 1.5, 2.0, 2.7, 5.5]

erg_bins = zip(erg_bins[:-1], erg_bins[1:])

# c1 = ROOT.TCanvas()
# c1.Divide(2, 2)
# c1.SetLogy()

total_Al_events = 0
total_DU_events = 0

histos = []
i=1


total_Al_histDP = TH1F(20, 180, binwidths=10)
total_DU_histDP = TH1F(20, 180, binwidths=10)

total_Al_histDP.Project(Al_treeDP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', weight=1.0/2/Al_pulses_DP_doubles)
total_DU_histDP.Project(DU_treeDP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', weight=1.0/2/DU_pulses_DP_doubles)

x_DU=total_DU_histDP.bincenters[0]
y_DU = total_DU_histDP.binvalues
y_DU_err = total_DU_histDP.binerrors

x_Al=total_Al_histDP.bincenters[0]
y_Al = total_Al_histDP.binvalues
y_Al_err = total_Al_histDP.binerrors


opts = {"capsize":3, "color":"black", "markersize":5, "elinewidth":1}

plt.errorbar(x_DU,y_DU, y_DU_err, linestyle='None', marker=".", label="DU target",**opts )
plt.errorbar(x_Al,y_Al, y_Al_err, linestyle='None', marker="x", label="Al target (noise estimate)", **opts)

plt.yscale('log')
plt.ylabel("counts per pulse")
plt.xlabel(r"$\theta_{nn}$ [deg]")
plt.legend(bbox_to_anchor=(0.4, 0.9, 1., .102),framealpha=1)
plt.subplots_adjust(bottom=0.13, top=0.8)
print sum(y_Al)/sum(y_DU)*100

plt.savefig('/Users/jeffreyburggraf/Pictures/Noise.png', bbox_inches='tight')


plt.show()



# total_Al_histDP.Draw(make_new_canvas=False)
# total_DU_histDP.Draw("same", make_new_canvas=False)
# total_Al_histDP.SetMaximum(max(total_DU_histDP.binvalues))
#
# total_Al_histDP.SetStats(0)
#
#
#
# mt2.thesis_plot([total_DU_histDP, total_Al_histDP], big_font=0.065)
#
# leg = ROOT.TLegend()
# leg.AddEntry(total_Al_histDP, 'Al target (noise estimate)', 'f')
# leg.AddEntry(total_DU_histDP, 'DU target')
# leg.Draw()

#
# for E1, E2 in erg_bins:
#     pad = c1.cd(i)
#
#     pad.SetGridy()
#
#     Al_histDP = TH1F(20, 180, binwidths=10)
#     DU_histDP = TH1F(20, 180, binwidths=10)
#
#     histos.extend([Al_histDP, DU_histDP])
#
#     cut = mt.cut_rangeAND([E1, E2], '0.5*(neutrons.coinc_hits[].erg[0] + neutrons.coinc_hits[].erg[1])')
#
#     Al_events = Al_histDP.Project(Al_treeDP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', cut=cut, weight=1.0/2/Al_pulses_DP_doubles)
#     DU_events = DU_histDP.Project(DU_treeDP_doubles, '180/3.1415*neutrons.coinc_hits[].coinc_theta', cut=cut, weight=1.0/2/DU_pulses_DP_doubles)
#     print "lllll"
#     scale = 0.75
#
#     Al_events=1
#
#     total_DU_events += DU_events
#     total_Al_events += Al_events
#
#     Al_histDP *= scale**2
#
#     total_Al_histDP += Al_histDP
#     total_DU_histDP += DU_histDP
#
#     ratio = (scale**2*Al_events/Al_pulses_DP_doubles)/(DU_events/DU_pulses_DP_doubles)
#
#     DU_histDP.Draw('hist E', make_new_canvas=False)
#     DU_histDP.SetStats(0)
#
#     DU_histDP.SetMinimum(0)
#     DU_histDP.SetLineWidth(3)
#
#     DU_histDP.SetLineColor(ROOT.kRed)
#     Al_histDP.Draw('hist same E')
#     Al_histDP.SetLineWidth(2)
#     Al_histDP.SetFillStyle(3114)
#     Al_histDP.SetFillColorAlpha(ROOT.kBlue,0.4)
#
#     DU_histDP.GetYaxis().SetTitle('nn_{uncorr}(#theta)')
#     DU_histDP.GetXaxis().SetTitle('#theta_{nn}')
#     DU_histDP.SetTitle('#frac{{#Sigma Al}}{{#Sigma DU}} = {0:.2f}, {1:.1f}<#bar{{E}}<{2:.1f} MeV'.format(ratio, E1, E2))
#     DU_histDP.SetTitleOffset(1.5)
#
#     mt2.thesis_plot([DU_histDP, Al_histDP], big_font=0.07)
#
#     DU_histDP.SetMinimum(0.5*min(Al_histDP.binvalues[np.where(Al_histDP.binvalues>0)]))
#     # DU_histDP.SetMinimum(1E-9)
#
#     pad.SetLogy(1)
#
#     i += 1
#     break
#
#
# c1.Draw()
#
#
# c2 = ROOT.TCanvas()
#
# # c2.SetGridy(1)
# c2.SetLogy(1)
# total_DU_histDP.Draw('hist E', make_new_canvas=False)
# total_Al_histDP.Draw('same hist E')
# total_DU_histDP.GetYaxis().SetTitle('nn_{uncorr}(#theta)')
# total_DU_histDP.GetXaxis().SetTitle('#theta_{nn} [deg]')
# total_DU_histDP.SetLineColor(ROOT.kRed)
# total_DU_histDP.SetLineWidth(3)
# total_DU_histDP.SetMinimum(1E-9)
# total_DU_histDP.SetStats(0)
#
# total_Al_histDP.SetLineWidth(2)
# total_Al_histDP.SetFillStyle(3114)
# total_Al_histDP.SetFillColorAlpha(ROOT.kBlue,0.4)
#

#
# mt2.thesis_plot([total_DU_histDP, total_Al_histDP], big_font=0.065)
#
#
#
# print('DU Events per pulse: {0:.2E}, Al events per pulse: {1:.2E}'.format(DU_events/DU_pulses_DP_doubles, scale**2*Al_events/Al_pulses_DP_doubles))
# print('Du pulses: {0} , Al pulses: {1}'.format(DU_pulses_DP_doubles, Al_pulses_DP_doubles))
# print('Total noise fraction: {:.3f}'.format(scale**2*total_Al_events/Al_pulses_DP_doubles/float(total_DU_events/DU_pulses_DP_doubles)))

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
                    if ___cmd___ == 'tb':
                        ___tb___ = ROOT.TBrowser()
                    else:
                        try:
                            exec (___cmd___, globals())
                        except:
                            print sys.exc_info()
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()


