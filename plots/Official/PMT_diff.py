import ROOT
import numpy as np
from TH1Wrapper import*
import mytools2 as mt2

import ROOT
import numpy as np
import mytools2 as mt2



dx = (2.54*30)/5

delta_t_mean = []
pos = []
errx = []
offset = 2.54*4


all_entries = []
positions = np.arange(offset,offset + dx*5,dx)
for p in positions:
    d1 = 30*2.54 - p
    d2 = p
    hist = TH1F(-10, 10, binwidths=0.25)
    print(p)

    entries =[]
    _x_ = d2 - 40.64
    sigma_t = np.e**(-np.abs(_x_/150))*2.5
    for delta in range(4000):
        _r_p = np.random.rand()
        _sig = np.random.randn()*sigma_t
        s1 = _sig*_r_p
        s2 = _sig*(1-_r_p)
        t1 = d1/(1.52E10)*1E9 + s1
        t2 = d2/(1.52E10)*1E9 + s2

        Dt = t2-t1
        entries.append(Dt)

    entries = np.array(entries) - 0.5

    pos.append(_x_)#- (2.54*30 + offset)/2.)
    delta_t_mean.append(np.mean(entries))

    all_entries.extend(entries)

    errx.append(sigma_t/np.sqrt(300))

pos = np.array(pos, dtype=np.float64)
delta_t_mean = np.array(delta_t_mean, dtype=np.float64)
errx = np.array(errx, dtype=np.float64)
erry = np.ones_like(errx)*0.5

gr = ROOT.TGraphErrors(len(pos), delta_t_mean, pos,errx,erry)

f = gr.Fit('pol1', 'S')


# projection_histo.Draw()

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt




font = {'family':'sans-serif',
        'sans-serif':['Helvetica'],
        'size': 25} # 18 for single figure, 30 for double

mpl.rc('font', **font)
mpl.rc('text', usetex=True)
mpl.rc("savefig", dpi=300)

project_hist = f.GetParams()[1]*np.array(all_entries) + f.GetParams()[0]
plt.figure(figsize=(10,10))
bins = np.arange(min(pos)-dx,max(pos)+dx, dx/15)
# bins += 0.5*(bins[1] - bins[0])
hist_y, hist_x = np.histogram(project_hist, bins=bins)
hist_x = hist_x[1:]
# plt.hist(project_hist, bins = np.arange(-40,40), fill=False,histtype='bar')
plt.errorbar(hist_x, hist_y,linewidth=1, drawstyle='steps-mid', elinewidth=1., mec='black', capsize=2, c='black')

plt.xticks(pos)
plt.grid()
plt.xlabel('Reconstructed distance from detector center [cm]')
plt.ylabel('Counts')

plt.savefig('/Users/jeffreyburggraf/Pictures/PMTDifference_hist.png', transparent=True)

plt.figure(figsize=(10,10))

plt.errorbar(delta_t_mean, pos, xerr=errx, yerr=erry, elinewidth=1.5, mec='black', capsize=4, c='black', linewidth=0, marker='d',
             label='measured data')

dx = (max(delta_t_mean)-min(delta_t_mean))/20.
fit_x = np.arange(min(delta_t_mean), max(delta_t_mean)+dx, dx)
fit_y = fit_x*f.GetParams()[1] + f.GetParams()[0]
plt.plot(fit_x, fit_y, c='black', ls='--', label='lin. fit: $y = ({p1:.1f} \pm {p1e:.1f})x+({p0:.0f} \pm {p0e:.1f})$'
         .format(p1=f.GetParams()[1], p1e=f.ParError(1), p0=f.GetParams()[0], p0e=f.ParError(0)))

# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xlabel(r'PMT timing difference $\pm$ SD [ns]')
plt.ylabel(r'Source distance from detector center [cm]')

plt.minorticks_on()
plt.yticks(np.arange(-40, 40+10, 10))
plt.xticks([-4,-2,0,2,4])

plt.grid()
plt.legend(loc='upper left', fontsize=25)
plt.ylim(-40, 55)
plt.subplots_adjust(left = 0.15)

plt.savefig('/Users/jeffreyburggraf/Pictures/PMTDifference.png', transparent=True)

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