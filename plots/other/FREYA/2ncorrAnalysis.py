import numpy
import ROOT
import mytools as mt
import re
from TH1Wrapper import TH1F

f = open('/Users/jeffreyburggraf/FREYA/MyFREYA/2nCorr/result.txt')

class Event:
    def __init__(self, erg, dir):
        self.erg = erg
        self.dir = dir

class Pulse:
    def __init__(self, photon_erg):
        self.photon_erg = [photon_erg]
        self.events = [[]]
        self.index = 0
        self.corr_theta = []
        self.uncorr_theta = []

    def new_pulse(self, photon_erg):
        self.photon_erg.append(photon_erg)
        self.events.append([])
        self.index += 1

    def add_event(self, data):
        erg = data[0]
        dir = data[1:]
        self.events[self.index].append(Event(erg, dir))

    def process(self):

        for event_list in self.events:
            for i in range(len(event_list) - 1):
                for j in range(i+1, len(event_list)):
                    e1 = event_list[i]
                    e2 = event_list[j]
                    self.corr_theta.append(180/3.1415*mt.vectorAngle(e1.dir, e2.dir))

        for i in range(len(self.events[0])):
            for j in range(len(self.events[1])):
                e1 = self.events[0][i]
                e2 = self.events[1][j]
                theta = 180/3.1415*mt.vectorAngle(e1.dir, e2.dir)
                if theta<2:
                    print(e1.erg, e1.dir), e2.erg, e2.dir
                self.uncorr_theta.append(theta)


line = ''
values = []
def get_line():
    global line,values
    line = f.readline()
    values = map(float, line.split())
    return line


def run(max_events = None):
    pulses = []
    pulse = None


    get_line()
    n = 0
    while True:
        if max_events:
            if n>max_events:
                return pulses

        if len(values) == 1:
            if pulse is None:
                pulse = Pulse(values[0])
            elif pulse.index == 1:
                pulse.process()
                n += len(pulse.uncorr_theta)
                pulses.append(pulse)
                pulse = Pulse(values[0])

            else:
                pulse.new_pulse(values[0])

        else:
            pulse.add_event(values)

        if not get_line():
            return pulses


pulses = run(40000)

histSP = TH1F(0,180,binwidths=5)
histDP = TH1F(0,180,binwidths=5)

for pulse in pulses:
    for theta in pulse.corr_theta:
        histSP.Fill(theta)

    for theta in pulse.uncorr_theta:
        histDP.Fill(theta)

histSP.normalize()
histDP.normalize()

histSP.SetMinimum(0)

histSP.Draw('hist')
histDP.Draw('same hist')


hist_W = histSP/histDP
hist_W.SetMinimum(0)
hist_W.Draw('hist')




if __name__ == "__main__":
    import ROOT as ROOT
    from multiprocessing import Process, Queue
    import time, sys, os


    def input_thread(q, stdin):
        while True:
            print ('ROOT: ')
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
                        print (sys.exc_info())
                time.sleep(0.01)
                ___g___()
        except(KeyboardInterrupt, SystemExit):
            ___input_p___.terminate()
