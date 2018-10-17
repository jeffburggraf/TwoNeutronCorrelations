import ROOT
import numpy as np
import mytools2 as mt2
import re
from TH1Wrapper import *

s="""(32400*ex2**2*(-((x2*(x1**2 + y1**2 + z1**2)*(x1*x2 + y1*y2 + z1*z2))/
     -           ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + 
     -        x1/Sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/
     -   (Pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/
     -        ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + 
     -  (32400*ex1**2*(-((x1*(x1*x2 + y1*y2 + z1*z2)*(x2**2 + y2**2 + z2**2))/
     -           ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + 
     -        x2/Sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/
     -   (Pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/
     -        ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + 
     -  (32400*ey2**2*(-((y2*(x1**2 + y1**2 + z1**2)*(x1*x2 + y1*y2 + z1*z2))/
     -           ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + 
     -        y1/Sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/
     -   (Pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/
     -        ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + 
     -  (32400*ey1**2*(-((y1*(x1*x2 + y1*y2 + z1*z2)*(x2**2 + y2**2 + z2**2))/
     -           ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + 
     -        y2/Sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/
     -   (Pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/
     -        ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + 
     -  (32400*ez2**2*(-(((x1**2 + y1**2 + z1**2)*z2*(x1*x2 + y1*y2 + z1*z2))/
     -           ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + 
     -        z1/Sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/
     -   (Pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/
     -        ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + 
     -  (32400*ez1**2*(-((z1*(x1*x2 + y1*y2 + z1*z2)*(x2**2 + y2**2 + z2**2))/
     -           ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + 
     -        z2/Sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/
     -   (Pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/
     -        ((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))"""

ROOT.gROOT.ProcessLine(".L twontree.h")

treeSP, _ = mt2.NCorrRun("SP","DU",generate_dictionary=0,Forward = True).neutrons_doubles_tree

def get_error(x1,y1,z1,x2,y2,z2, ex1,ey1,ez1, ex2,ey2,ez2):
    result =  (32400*ex2**2*(-((x2*(x1**2 + y1**2 + z1**2)*(x1*x2 + y1*y2 + z1*z2))/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + x1/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/(np.pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + (32400*ex1**2*(-((x1*(x1*x2 + y1*y2 + z1*z2)*(x2**2 + y2**2 + z2**2))/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + x2/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/(np.pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + (32400*ey2**2*(-((y2*(x1**2 + y1**2 + z1**2)*(x1*x2 + y1*y2 + z1*z2))/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + y1/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/(np.pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + (32400*ey1**2*(-((y1*(x1*x2 + y1*y2 + z1*z2)*(x2**2 + y2**2 + z2**2))/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + y2/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/(np.pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + (32400*ez2**2*(-(((x1**2 + y1**2 + z1**2)*z2*(x1*x2 + y1*y2 + z1*z2))/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + z1/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/(np.pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))) + (32400*ez1**2*(-((z1*(x1*x2 + y1*y2 + z1*z2)*(x2**2 + y2**2 + z2**2))/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))**1.5) + z2/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2)))**2)/(np.pi**2*(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))

    return np.sqrt(result)


b_wdth = 3
err_hist = TH1F(20,180, binwidths=b_wdth )
n_fills = TH1F(20,180, binwidths=b_wdth )

errs = []
for evt in treeSP:
    hit = evt.neutrons.coinc_hits[0]
    x1,y1,z1 = hit.x[0], hit.y[0], hit.z[0]
    x2,y2,z2 = hit.x[1], hit.y[1], hit.z[1]

    det1 = hit.det[0]
    det2 = hit.det[1]

    ex1 = 2.54*3*np.sin((np.pi*det1)/180.)
    ey1 = 2.54*3*np.cos((np.pi*det1)/180.)

    if det1 in [30,330]:
        ez1 = 2.54*10
    else:
        ez1 = 2.54 *5
    if det2 in [30,330]:
        ez2 = 2.54*10
    else:
        ez2 = 2.54 *5

    ex2 = 2.54 * 3 * np.sin((np.pi * det2) / 180.)
    ey2 = 2.54 * 3 * np.cos((np.pi * det2) / 180.)

    err = (get_error(x1,y1,z1,x2,y2,z2,ex1,ey1,ez1, ex2,ey2,ez2))

    theta = (180*np.arccos((x1*x2 + y1*y2 + z1*z2)/np.sqrt((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/np.pi
    print theta,err
    err_hist.Fill(theta,err)
    n_fills.Fill(theta)

err_hist /= n_fills

print (list(err_hist.binvalues))
print (list(err_hist.__binLeftEdges__))

err_hist.Draw("hist p")
err_hist.SetStats(0)
mt2.thesis_plot([err_hist],big_font=0.045)

err_hist.GetXaxis().SetTitle("Reconstructed opening angle [degrees]")
err_hist.GetYaxis().SetTitle("Mean uncertainty in opening angle [degrees]")
err_hist.SetMarkerStyle(22)
n_fills.Draw()


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


