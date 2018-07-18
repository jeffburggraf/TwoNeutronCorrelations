import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2

import ROOT
import numpy as np
import mytools2 as mt2



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
    y.append(np.mean(entries)-0.5)
    erry.append(np.std(entries)/np.sqrt(len(entries)))

x = np.array(x, dtype=np.float64)
y = np.array(y, dtype=np.float64)
erry = np.array(erry, dtype=np.float64)

gr = ROOT.TGraphErrors(len(x), x, y,np.ones_like(erry)*0.5, erry)

f = gr.Fit('pol1', 'S')


import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

plt.figure(figsize=(10,10))


font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size': 30} # 18 for single figure, 30 for double

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc("savefig", dpi=300)

plt.errorbar(x, y, yerr=erry, elinewidth=1.5, mec='black', capsize=4, c='black', linewidth=0, marker='d',
             label='measurement')

dx = (x[-1]-x[0])/20.
fit_x = np.arange(x[0], x[-1]+dx, dx)
fit_y = fit_x*f.GetParams()[1] + f.GetParams()[0]
plt.plot(fit_x, fit_y, c='black', ls='--', label='lin. fit: $y = (0.13 \pm 0.02)x+(0.16 \pm 0.5)$')

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'Distance from detector center [cm]')
plt.ylabel(r'PMT timing difference [ns]')

plt.minorticks_on()
plt.xticks(np.arange(-40, 40+10, 10))
plt.yticks([-4,-2,0,2,4])

plt.grid()
plt.legend(loc='upper left', fontsize=25)
plt.ylim(min(y)*1.15, 6.75)

# plt.savefig('/Users/jeffreyburggraf/PycharmProjects/2nCorrPhysRev/PMTDifference.png', transparent=True)

# plt.savefig('PMTDifference.png', transparent=True)
plt.show()


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