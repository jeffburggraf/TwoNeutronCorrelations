import os,time,sys
import numpy as np
import warnings
import mytools as mt
from numbers import Number
import TH1Wrapper

import re

try:
    import ROOT
except:
    warnings.warn("ROOT module not loaded. Not all features are available.")

class NCorrRun:
    files = []
    def __init__(self,SP_or_DP, name_target, generate_dictionary = False , Forward = False):
        self.generate_dictionary = generate_dictionary

        if (sys.platform=="darwin"):
            if not Forward:
                base_directory = "/Volumes/JeffMacEx/2nCorrData/Production/"
            else:
                base_directory = "/Volumes/JeffMacEx/2nCorrData/Production/Forward/"
        else:
            assert False, "Need to fix this"

        assert isinstance(SP_or_DP,str)
        assert "SP" in SP_or_DP or "DP" in SP_or_DP
        assert isinstance(name_target, str)
        assert os.path.isdir(base_directory) , "Plug in your external hard drive!"

        self.file_dir = "{0}{1}_{2}/".format(base_directory, SP_or_DP ,name_target)

        self.file_name_prefix ="{0}_{1}".format( SP_or_DP, name_target)

        assert os.path.isdir(self.file_dir ), "could not find dir {}".format(self.file_dir )

        self.ROOT_file_names_dict = {}

    def __getfile__(self, root_file_name):

        if self.generate_dictionary and len(self.ROOT_file_names_dict) == 0:
            os.getenv('DYLD_LIBRARY_PATH')
            ROOT.gSystem.CompileMacro("/Users/jeffreyburggraf/root-6.12.04/macros/twontree.h")

            # ROOT.gROOT.ProcessLine(".L twontree.h+")

        root_file_name = self.file_name_prefix + root_file_name
        if root_file_name not in self.ROOT_file_names_dict:
            abs_path = self.file_dir + root_file_name+".root"
            assert os.path.exists(abs_path), "cannot find file: {}".format(abs_path)
            self.ROOT_file_names_dict[root_file_name]=ROOT.TFile(abs_path)

        return self.ROOT_file_names_dict[root_file_name]

    def n__get_tree_and_puleses__(self,root_file_name):

        file = self.__getfile__(root_file_name)

        NCorrRun.files.append(file)
        for tkey in file.GetListOfKeys():
            if tkey.GetName()=="tree":break
        else:
            assert False, "Couldn't find key named 'tree' in file. The keys available are: {}".format([k for k in  file.GetListOfKeys()])

        tree = file.Get("tree")

        assert hasattr(tree, "PulseNumber"), "No PulseNumber branch in tree. This is needed for finding the number of pulses looked at!"

        tree.GetEntry(tree.GetEntries()-1)

        n_pulses = tree.PulseNumber
        tree.GetEntry(0)
        return tree, float(n_pulses)

    @property
    def neutrons_doubles_tree(self):
        return self.n__get_tree_and_puleses__("_neutrons_coinc")
    @property
    def neutrons_singles_tree(self):
        return self.n__get_tree_and_puleses__("_neutrons_singles")

    @property
    def photons_doubles_tree(self):
        return self.n__get_tree_and_puleses__("_photons_coinc")

    @property
    def photons_singles_tree(self):
        return self.n__get_tree_and_puleses__("_photons_singles")
    
    @property
    def noise_doubles_tree(self):
        return self.n__get_tree_and_puleses__("_noise_coinc")

    @property
    def noise_singles_tree(self):
        return self.n__get_tree_and_puleses__("_noise_singles")

    @property
    def all_doubles_tree(self):
        return self.n__get_tree_and_puleses__("_all_coinc")

    @property
    def all_singles_tree(self):
        return self.n__get_tree_and_puleses__("_all_singles")

    def get_nPulses_for_run(self, runnum):
        tree,_ = self.neutrons_doubles_tree
        run_start = tree.Project("","","RunNumber<{0}".format(runnum))
        run_end = run_start + tree.Project("","","RunNumber=={0}".format(runnum),"",tree.GetEntries(), run_start) - 1

        tree.GetEntry(run_start)
        pulse_start = tree.PulseNumber
        tree.GetEntry(run_end)
        pulse_end = tree.PulseNumber
        return pulse_end-pulse_start



def get_rate(tree,cut="",ndigits=3):
    assert isinstance(tree,ROOT.TTree)

    assert hasattr(tree,"PulseNumber")
    tree.GetEntry(tree.GetEntries() - 1)
    n_pulses = float(tree.PulseNumber)
    tree.GetEntry(0)

    n_events = tree.Project("","",cut)
    err = np.sqrt(n_events)/n_pulses
    rate = n_events/n_pulses

    string = "{0} +/- {1}".format(mt.round(rate,ndigits),mt.round(err,1))
    return string, rate,err


def get_weighted_avg_err(values,weights,errors):
    subterm = sum(values*weights)/sum(weights)**2
    varryterms = values/sum(weights)
    partialsums = [(vt-subterm)**2*sigma**2 for vt,sigma in zip(varryterms,errors)]
    return np.sqrt(sum(partialsums))

def get_divided_Error(Num,Denom,Numerror,DenomError):
    return np.sqrt(
        (DenomError ** 2 * Num ** 2 + Numerror ** 2 * Denom ** 2) / (
                Denom ** 4))

def leg_get_leg_position(leg):
    ROOT.gPad.Update()
    ROOT.gPad.Modified()
    print  map(lambda x:round(x,2),[leg.GetX1(), leg.GetY1(),leg.GetX2(),leg.GetY2()])

def set_leg_position(x1,y1,x2,y2,leg):
    ROOT.gPad.Update()
    ROOT.gPad.Modified()
    leg.SetX1(x1), leg.SetY1(y1), leg.SetX2(x2), leg.SetY2(y2)

def thesis_plot(*plots):
    if not hasattr(plots,"__iter__"):
        plots = [plots]

    for graph in plots:
        graph.GetXaxis().CenterTitle()
        graph.GetYaxis().CenterTitle()

        graph.GetXaxis().SetTitleOffset(1.5)
        graph.GetYaxis().SetTitleOffset(1.5)

        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.13)

import warnings

class row:
    def __init__(self, col_spans , label = None, rowspan = 1):
        self.width = sum(col_spans)

        self.entries = [""]*self.width

        self.ientry = 0

        self.label = label
        self.rowspan =rowspan

        self.is_full = False

    def add_entry(self, value):
        assert self.ientry is not None, "Cannot insert values sequentially after table.set_value() has been used!"
        assert self.ientry < self.width, "Too many entries added to row. Max = {0}".format(self.width)
        self.entries[self.ientry] = "{0}".format(value)

        self.ientry += 1

        if self.ientry == self.width:
            self.is_full = True

    def __get_text__(self):
        s = ""
        if self.label:
            s += "! {0}\n".format(self.label)

        s += "| " + " || ".join(self.entries)

        return s

class table:
    def __init__(self):
        self.column_labels = []
        self.row_labels = []

        self.colomn_spans = []

        self.rows = []

        self.current_row = None

    def create_colomn(self, label, span = 1):
        assert isinstance(span, int)
        assert self.current_row is None, "Cannot create new columns after creating a row!"

        self.column_labels.append("{}".format(label))
        self.colomn_spans.append(span)

    def new_row(self, label = None, row_span = 1):
        assert self.column_labels, "Add some columns before adding a row!"

        new_row = row(self.colomn_spans, label, row_span)

        if self.current_row is not None:
            if not self.current_row.is_full:
                if not hasattr(self,"stop_row_not_full_warning"):
                    self.stop_row_not_full_warning  = True
                    warnings.warn("Beginning new row before current row is completely filled! Further warnings suppressed. ")

                self.current_row.add_entry("")

        if label is not None:
            assert len(self.row_labels)== len(self.rows), "If using row labels, all rows should have a label!"
            self.row_labels.append(label)

        self.current_row = new_row
        self.rows.append(self.current_row)

        return self.current_row

    def __get_column_label_lines__(self):
        column_labels = []
        for label,span in zip(self.column_labels,self.colomn_spans):
            if span==1:
                column_labels.append(label)
            else:
                column_labels.append("colspan = {0} | {1}".format(span,label))
        s = ""
        if self.row_labels:
            s += "! !"

        s += "! " + " !! ".join([s for s in column_labels])
        return s

    def set_value(self,column_label,row_label, value):
        column_label = str(column_label)
        row_label = str(row_label)

        assert len(self.colomn_spans) == sum(self.colomn_spans), "All columns must have a width of exactly one to use this method. "
        assert self.row_labels, "Must using row labels to use this method!"
        assert column_label in self.column_labels, "Cannot find column with label '{0}'".format(column_label)
        assert row_label in self.row_labels, "Cannot find row with label '{0}'".format(row_label)

        for row in self.rows:
            if row.label == row_label:
                break
        else:
            assert False, "Invalid row label"

        row.ientry = None
        row.entries[self.column_labels.index(column_label)] = "{0}".format(value)

    def __repr__(self):
        lines = [self.__get_column_label_lines__()]

        for row in self.rows:
            lines.append(row.__get_text__())

        result = "\n|-\n".join(lines)
        result ='{| class="wikitable"\n'  + result + "\n|}"

        return result


runs = {"DU": [6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6542, 6543, 6544,
               6545, 6546, 6547, 6548, 6549, 6553, 6554, 6556, 6557, 6558, 6559,
               6560, 6561, 6568, 6570, 6571, 6572, 6573, 6587, 6588, 6591, 6592,
               6600, 6601, 6603, 6608, 6609, 6610, 6611, 6612, 6613, 6614, 6615,
               6616, 6617, 6625, 6626, 6627, 6628, 6629, 6630, 6634, 6635, 6636,
               6640],
        "D2O": [6541, 6551, 6565, 6586, 6607, 6621, 6642],
        "Th": [6622, 6623, 6624, 6638, 6639],
        "Al": [6604, 6552, 6569, 6598, 6620, 6632, 6643]
        }

runs_good = {
    "DU": [6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6542, 6543, 6544, 6545, 6546, 6547, 6548, 6549, 6553, 6554,
           6556, 6557, 6558, 6568, 6570, 6571, 6572, 6573, 6591, 6592, 6600, 6601, 6603, 6608, 6609, 6610, 6611, 6612,
           6613, 6615, 6616, 6617, 6627, 6634, 6635, 6636, 6640],
    "D2O": [6541, 6565, 6607, 6621, 6642],
    "Th": [6622, 6638, 6639]
        }

runs_bad ={"DU":[6559, 6560, 6561, 6587, 6588, 6614, 6625, 6626, 6628, 6629, 663],
           "D2O":[6551, 6586],
           "Th":[6623, 6624]
           }


class draw_all():
    def __init__(self):
        self._Drawables = [[]]
        self._options = [[]]
        self._funcs = [[]]

    def add_draw_ables(self, drawable, option = "", func = lambda x:x):

        assert hasattr(drawable, "Draw"), "drawables argument must have Draw method!"

        self._Drawables[-1].append(drawable)
        self._options[-1].append(option)
        assert hasattr(func, "__call__")

        self._funcs[-1].append(func)

    def end_plot(self):
        self._Drawables.append([])
        self._options.append([])
        self._funcs.append([])

    def Draw(self):
        
        self._Drawables = filter(lambda x :len(x)!=0, self._Drawables)
        self._options = filter(lambda x :len(x)!=0, self._options)

        num_things_total = len(self._Drawables)


        assert not hasattr(self, "self.canvii"), "Cant draw twice!"

        self.canvii = []
        pad_index = overflow_pad_index = None

        num_things_left = num_things_total
        num_things_drawn = 0

        for draw_ables,options,funcs in zip(self._Drawables,self._options, self._funcs):
            num_things_left = num_things_total - num_things_drawn

            if pad_index==None or pad_index==overflow_pad_index:

                self.canvii.append(ROOT.TCanvas())
                print "making new TCanvas()"
                print "num_things_left: {0}".format(num_things_left)


                c = self.canvii[-1]

                if num_things_left==1:
                    overflow_pad_index = 2
                elif num_things_left == 2:
                    c.Divide(2)
                    overflow_pad_index = 3
                elif 2<num_things_left<=4:
                    overflow_pad_index = 5
                    c.Divide(2,2)
                elif 5<=num_things_left<=6:
                    c.Divide(3,2)
                    overflow_pad_index = 7
                elif 7<=num_things_left<=9:
                    c.Divide(3,3)
                    overflow_pad_index = 10
                else:
                    overflow_pad_index = 13
                    c.Divide(4,3)

                pad_index = 1
                print "overflow_pad_index: {0}\n".format(overflow_pad_index)
                print "pad_index: {0}\n".format(pad_index)


            c.cd(pad_index)
            assert enumerate(zip(draw_ables, options))

            for index, (drw,opt,func) in enumerate(zip(draw_ables, options,funcs)):
                kwargs = {}
                if isinstance(drw,TH1Wrapper.hist):
                    kwargs["make_new_canvas"] = False

                if index==0:
                    drw.Draw(opt, **kwargs)
                else:
                    drw.Draw(opt + "same", **kwargs)

                func(drw)


            pad_index += 1
            num_things_drawn +=1

def median(hist, n):
    assert isinstance(hist, TH1Wrapper.TH1F)
    assert isinstance(n,int)

    value = sum(hist.binvalues)/ float(n)

    sum_ = 0

    results = [hist.__binLeftEdges__[0][0]]
    for index,bin in enumerate(hist.binvalues):

        sum_ += bin

        if sum_ >= value:
            results.append(hist.__binLeftEdges__[0][index])
            sum_ = 0

    return results + [hist.__binRightEdges__[0][-1]]




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





