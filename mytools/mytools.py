import os,time,sys
import numpy as np
from itertools import izip_longest
from scipy import ndimage
from scipy import signal
import warnings,array,re
from collections import OrderedDict
from shutil import copyfile
from operator import itemgetter
from itertools import combinations
import numbers
from decimal import Decimal
from scipy.special import legendre

""" TTop"""
try:

    import ROOT
    pass
except:
    warnings.warn("ROOT module not loaded. Not all features are available.")
    pass

def Carr(list):
    return np.array(list,dtype=np.float32)

class rebin():
    def __init__(self,TH1F,newBinLeftEdges,AllowZeroBin=False):
        hlist  = HistToList(TH1F)

        assert sorted(newBinLeftEdges)==newBinLeftEdges,"Bins must be in increasing order!"

        self.binLeftEdges = newBinLeftEdges
        self.goodBins=np.arange(0,len(self.binLeftEdges)-1)
        self.HasBeenNormalized=False

        newBinIndex = 0

        _Xvalues_ = []
        _XWeights_=[]
        _YErrs_ = []

        self.x=np.zeros(len(newBinLeftEdges)-1,dtype=np.float32)
        self.y=np.zeros(len(newBinLeftEdges)-1,dtype=np.float32)
        self.erry = np.zeros(len(newBinLeftEdges)-1,dtype=np.float32)
        self.errx = np.zeros(len(newBinLeftEdges)-1,dtype=np.float32)
        self.bin_widths = np.array([i2-i1 for i1,i2 in zip(newBinLeftEdges[:-1],newBinLeftEdges[1:])])

        for i,(binCenter,value) in enumerate(hlist.xyTuple):
            if binCenter<newBinLeftEdges[0]:continue

            if binCenter >= newBinLeftEdges[newBinIndex+1]:
                if sum(_XWeights_)!=0:

                    _Xvalues_ = np.array(_Xvalues_)
                    _XWeights_ = np.array(_XWeights_)
                    _YErrs_ = np.array(_YErrs_)
                    x_mean_ = np.average(_Xvalues_,weights=_XWeights_)


                    self.x[newBinIndex] = x_mean_
                    self.y[newBinIndex] = sum(_XWeights_)  #/len(filter(lambda x:x>0,_XWeights_))
                    self.errx[newBinIndex] = np.sqrt(np.average((_Xvalues_-x_mean_)**2,weights=abs(_XWeights_)))
                    self.erry[newBinIndex] = np.sqrt(sum(_YErrs_*_YErrs_))

                else:
                    self.x[newBinIndex]= np.average(_Xvalues_)
                    self.y[newBinIndex]=0
                    self.errx[newBinIndex]=0
                    self.erry[newBinIndex]=0


                if (newBinIndex+1)==len(newBinLeftEdges) or binCenter>newBinLeftEdges[-1]:
                    break

                _Xvalues_ = []
                _XWeights_ = []
                _YErrs_ = []
                newBinIndex += 1

            _Xvalues_.append(binCenter)
            _XWeights_.append(value)
            _YErrs_.append(hlist.erry[i])

        else:
            if sum(_XWeights_) == 0:
                self.x[newBinIndex] = np.average(_Xvalues_)
        self.__set_max__()

    def __len__(self):
        return len(self.goodBins)

    def remove_bins_above_rel_err(self,rel_err = 0.3):
        self.__select__(np.where(self.y!=0))
        self.__select__(np.where(np.abs(self.erry/self.y)<rel_err))

    def multiplyByxyvalues(self,x,y):
        binminmax = zip(self.binLeftEdges[0:-1],self.binLeftEdges[1:])
        binminmax = [binminmax[i] for i in self.goodBins]
        X,Y, = list(x),list(y)
        for ixy, (_x_,_y_) in enumerate(zip(X,Y)):
            for ibin, (_min, _max) in enumerate(binminmax):
                if _min<=_x_<_max:
                    self.y[ibin] *= _y_
                    self.erry[ibin] *= _y_
                    break
        self.__set_max__()


    def reduce_to_same_bins(self,other):
        assert isinstance(other,rebin)

        if len(self) < len(other):
            smaller_set = self.goodBins
            larger_set = other.goodBins
            to_reduce= other

        elif len(self)>len(other):
            smaller_set = other.goodBins
            larger_set = self.goodBins
            to_reduce = self

        else:
            return 0

        good_bins = np.where([(i in smaller_set) for i in larger_set])

        to_reduce.__select__(good_bins)

        return len(other)-len(self)


    def __set_max__(self):
        self.max = max(self.y)

    def getTGraph(self):
        g = ROOT.TGraphErrors(len(self),self.x,self.y,self.errx,self.erry)
        return g


    def normalize(self,normalize_for_binWidth=True):
        if normalize_for_binWidth:

            self.y /= self.bin_widths  # Scale for each binwidth
            self.erry/=self.bin_widths

            norm = sum(self.y*self.bin_widths) # norm is effectively the area of the function integrated over all bins.
            self.y /= norm
            self.erry /=norm
            self.HasBeenNormalized = True
            self.__set_max__()
            return self

        norm = sum(self.y*self.bin_widths) # norm is effectively the area of the function integrated over all bins.
        self.y/=norm
        self.erry/=norm

        self.HasBeenNormalized = True
        self.__set_max__()


    def __select__(self,indicies):
        self.x = self.x[indicies]
        self.y = self.y[indicies]
        self.errx = self.errx[indicies]
        self.erry = self.erry[indicies]
        self.bin_widths = self.bin_widths[indicies]
        self.goodBins = self.goodBins[indicies]
        self.__set_max__()

        if self.HasBeenNormalized==True:
            self.normalize(False)


    def Divide(self,other):
        if isinstance(other,rebin):
            assert self.binLeftEdges==other.binLeftEdges

            assert not all(other.y==0), "All bins in the denominator are zero!"

            nbins_change = self.reduce_to_same_bins(other)
            if (nbins_change!=0):
                warnings.warn("Added/removed {0} bins from rebin obj: {1}".format(nbins_change, self))

            non_zero_bins = np.where(other.y>0)
            self.__select__(non_zero_bins)
            other.__select__(non_zero_bins)

            self.x = (self.x*self.y+other.x*other.y) / (self.y+other.y)
            self.errx = np.sqrt(self.errx**2+other.errx**2)


            self.erry = np.sqrt((other.erry**2*self.y**2 + self.erry**2*other.y**2)/(other.y**4))
            self.y = self.y / other.y
            self.__set_max__()
            return self

        elif isinstance(other,np.ndarray):
            assert len(other)==len(self)

            non_zero_bins = np.where(other>0)
            self.__select__(non_zero_bins)
            other  = other[non_zero_bins]

            self.y/=other
            self.erry/=other

            self.__set_max__()

            return self

        else:
            raise NotImplementedError


    def Add(self,other,c1=1,constErr=None):
        if isinstance(other,rebin):
            assert self.binLeftEdges==other.binLeftEdges

            chang_in_bins = self.reduce_to_same_bins(other)
            if chang_in_bins!=0:
                warnings.warn("Removed {1} bins from rebin obj: {0}".format(self,chang_in_bins))

            if  any((self.y + other.y)==0):
                for i in range(len(self.x)):
                    if self.y[i]==0 and other.y[i]==0:
                        self.x[i] = 0.5*(self.x[i] + other.x[i])
                    else:
                        self.x[i]= (self.x[i]*self.y[i] + other.x[i]*other.y[i]) / (self.y[i] + other.y[i])
            else:
                self.x = (self.x*self.y + other.x*other.y) / (self.y + other.y)

            self.y += c1*other.y

            self.errx = np.sqrt(other.errx**2+self.errx**2)
            self.erry = np.sqrt((c1*other.erry)**2 + self.erry ** 2)

            self.__set_max__()
            return self

        elif isinstance(other, numbers.Number):
            self.y += other

            if constErr:
                self.erry = np.sqrt(constErr ** 2 + self.erry ** 2)

            self.__set_max__()

            return self

        elif isinstance(other,np.ndarray) and len(self)==len(other):
            self.y*=other
            self.erry*=other
            self.__set_max__()
            return self

        else:
            raise NotImplementedError

    def __iter__(self):
        return iter(self.y)

    def Multiply(self,other,constErr=None):
        if isinstance(other, rebin):
            assert self.binLeftEdges == other.binLeftEdges

            if  (self.y + other.y)==0:
                self.x = 0.5*(self.x + other.x)
            else:
                self.x = (self.x * self.y + other.x * other.y) / (self.y + other.y)

            self.y = self.y*other.y
            self.errx = np.sqrt(other.errx ** 2 + self.errx ** 2)
            self.erry = np.sqrt((self.erry*other.y)**2+(self.y*other.erry)**2)
            self.__set_max__()
            return self

        elif isinstance(other, numbers.Number):
            self.y *= other
            if constErr:
                self.erry = np.sqrt((self.erry * other) ** 2 + (self.y * constErr) ** 2)
            else:
                self.erry *= other
            self.__set_max__()
            return self
        else:
            raise NotImplementedError

    def __div__(self, other):
        return self.Divide(other)

    def __rdiv__(self, other):
        assert isinstance(other,numbers.Number)
        other = float(other)
        self.y = other/self.y
        self.erry = other*self.erry/self.y**2
        self.__set_max__()
        return self

    def __add__(self, other):
        return self.Add(other)

    def __sub__(self, other):
        return self.Add(other,-1)

    def __neg__(self):
        self.y*=-1
        self.__set_max__()
        return self

    def __mul__(self, other):
        return self.Multiply(other)

def add_commas(num):
    s=str(num)
    s_rev=s[::-1]
    set_of_three=re.findall("[0-9]{3}",s_rev)

    remaining=s_rev[3*(len(set_of_three)):][::-1]

    set_of_three = [_[::-1] for _ in set_of_three[::-1]]

    result =",".join(set_of_three)
    if remaining:
        result=remaining+","+result
    return result

angles=[30,54,78,102,126,150,210,234,258,282,306,330]  # most recent

angle2TDC1190Index = {angles[i]: (2*i+1,2*i+2) for i in range(len(angles))}

data_dir = "/Volumes/JeffMacEx/2nCorrData/Aug_2017/"


def TH1F(title, _min,_max,binwidth):
    return ROOT.TH1F(str(title),str(title),int(float(_max-_min)/binwidth),_min,_max)

class run_groups:
    D20 = [6541, 6551, 6565, 6586, 6607, 6621,6642]
    DU = [6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6542, 6543, 6544,
              6545, 6546, 6547, 6548, 6549, 6553, 6554, 6556, 6557, 6558, 6559,
              6560, 6561, 6568, 6570, 6571, 6572, 6573, 6587, 6588, 6591, 6592,
              6600, 6601, 6603, 6608, 6609, 6610, 6611, 6612, 6613, 6614, 6615,
              6616, 6617, 6625, 6626, 6627, 6628, 6629, 6630, 6634, 6635, 6636,
              6640]
    Th = [6622, 6623, 6624, 6638, 6639]
    Al= [6604,6552,6569,6598, 6620,6632,6643]
    Cosmics = [6562,6585,6618,6631,6641]
    All_runs = range(6531,6658)

def cartesian_product(*arrays):
    la = len(arrays)
    arr = np.empty([len(a) for a in arrays] + [la])
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, la)

def findsubsets(S,m):
    return set(combinations(S, m))

tree_dir="/Users/jeffburggraf/Trees/"
def get_run(run_number,newtree=True):
    dir="/Users/jeffburggraf/Trees/"
    if newtree:
        path="{0}newr{1}.root".format(dir,run_number)
    else:
        path = "{0}r{1}.root".format(dir, run_number)

    assert os.path.isfile(path)

    return ROOT.TFile(path)



def round(x,sig_figs=2):
    sig_figs=sig_figs-1
    if sig_figs<0:sig_figs=0
    a='%.{}E'.format(sig_figs) % Decimal(str(x))

    return a

def NCorr_ROOTFiles_list(top_directory="/Volumes/JeffMacEx/2nCorrData/"):
    assert os.path.isdir(top_directory), "Cannot find {0}. Try plugging in your external hard drive.".format(top_directory)
    if top_directory[-1]!="/":
        top_directory+="/"

    dirs=filter(os.path.isdir,[top_directory+filename for filename in os.listdir(top_directory)])
    files=[]

    for dir in dirs:
        for filename in os.listdir(dir):
            if ".root" in filename:
                files.append(dir+"/"+filename)

    return files

class MCNPTTreePrint():

    def __init__(self,tree):
        assert isinstance(tree,ROOT.TTree), "argument mut be an instance of TTree. "
        self.tree=tree
        self.line=[]
        self.attributes=["nps","event","bnknum","par","cell","issource","erg","deposited","isdetection","x","y","z","weight","term",
                         "mevee","mt","zaid","dirx","diry","dirz","time","incident_energy"]

    def Print(self,thingstoprint=["auto"],verbosity=1,Round=2):
        warns=0
        max_warns=4
        attributes=[]
        if verbosity==1:
            verbose=["nps","event","bnknum","par","cell","issource","erg","deposited","isdetection"]
        else:
            verbose = self.attributes

        for branch in self.attributes:
            if branch in verbose:
                attributes.append(branch)
        cout=[]
        if thingstoprint!=["auto"]:
            attributes+=thingstoprint



        for attr in attributes:
            if hasattr(self.tree,attr):
                value=getattr(self.tree,attr)
                if value // 1==value:
                    value=int(value)
                else:
                    value=round(value,Round)
                cout.append("{0}: {1}".format(attr,value))
            else:
                warns+=1
                assert warns <= max_warns, "The tree provided to class has too many missing branches. Is the tree from Polimiparse.c? "
                warnings.warn("Branch {} is not a branch of the provided tree")
        print "  ".join(cout)



def find_ROOT_TTree_from_TFile(tfile):
    assert(isinstance(tfile,ROOT.TFile))
    _file=tfile
    i = 0
    tree=None
    while bool(_file.GetListOfKeys().At(i)):
        try:
            _key = _file.GetListOfKeys().At(i).GetName()
        except:
            assert False, "Auto tree location functoin failed. Provide TTree name using the tree_name argument.  "

        if isinstance(_file.Get(_key), ROOT.TTree):
            tree = _file.Get(_key)
            break
        i += 1

    assert isinstance(tree,ROOT.TTree), "Can't find tree in root file \n {}".format(_file.ls())

    return tree









file_refs=[]
_quiet_NCorrTree_=0
def NCorrTree(runnumber,generate_dictionary = False, rundir="/Volumes/JeffMacEx/2nCorrData/",TTreeBaseName="2ncorr",tree_name="auto"):
    warnings.warn("Use mytools2 instead! this function is depreciated!")
    global _quiet_NCorrTree_
    if rundir[-1]!="/":rundir+="/"
    assert os.path.isdir(rundir), "Cant find simulations directory {}. Maybe plug in 'JeffMacEx'?".format(rundir)

    if isinstance(runnumber,str):
        file_name1 = "/Volumes/JeffMacEx/2nCorrData/Aug_2017/"+runnumber+".root"
        filename2 =  "/Volumes/JeffMacEx/2nCorrData/Aug_2017/Production/"+runnumber+".root"
        __dir__ = "/Volumes/JeffMacEx/2nCorrData/Aug_2017/ProductionTTrees/"

        if os.path.exists(filename2):
            file_name = filename2
        else:
            file_name = file_name1
            warnings.warn("Not using production trees. ")

        if not os.path.exists(file_name):
            print "\nOptions:\n".format(os.listdir(__dir__))

            assert False

        if generate_dictionary:
            if len(file_refs)==0:
                ROOT.gROOT.ProcessLine(".L twontree.h+")

        file_refs.append(ROOT.TFile(file_name))
        tree = file_refs[-1].Get("tree")
        tree.GetEntry(tree.GetEntries()-1)

        return tree, tree.PulseNumber

    directories=filter(lambda x:os.path.isdir(rundir+x),os.listdir(rundir))
    located_N_Matching_ROOT_Files=0
    found_directories=[]

    base_file_name=TTreeBaseName+str(runnumber)+".root"

    for dir in directories:
        dir=rundir+dir+"/"

        abs_path=dir+base_file_name
        if os.path.isfile(abs_path):
            _file= ROOT.TFile(abs_path)
            file_refs.append(_file)
            located_N_Matching_ROOT_Files+=1
            found_directories.append(dir)

    if not located_N_Matching_ROOT_Files:
        _quiet_NCorrTree_+=1
        if _quiet_NCorrTree_<=3:
            print "Cant find root tree for runnum {0}{1}.root".format(TTreeBaseName,runnumber)
        if _quiet_NCorrTree_==3:
            print "suppressing 'cant find' warnings.  "
        return None
    assert located_N_Matching_ROOT_Files<=1, "Duplicate ROOT files of name {0} found in: \n {1} ".format(base_file_name,found_directories)
    # assert located_N_Matching_ROOT_Files!=0, "Could not find root file named {}".format(base_file_name)

    if tree_name!="auto":
        tree=_file.Get(tree_name)
    else:
        tree=find_ROOT_TTree_from_TFile(_file)

    if isinstance(tree,ROOT.TTree):
        return tree
    else:
        return None


def NCorrMultiTree(run_numbers_list,TTreeBaseName="2ncorr",max_pulses=False,max_events_dict={},rundir="/Volumes/JeffMacEx/2nCorrData/",tree_name=None):
    '''
    :param run_numbers_list:
    :param max_events_list: list of max number of events to precess in each tree, respectively, e.g. [1234,234500] for runnumlist [6133,6144]
    :param rundir:
    :param TTreeBaseName:
    :param tree_name:
    :return:
    '''
    if rundir[-1]!="/":rundir+="/"
    assert os.path.isdir(rundir), "Cant find simulations directory. Maybe plug in 'JeffMacEx'?"

    directories=filter(lambda x:os.path.isdir(rundir+x),os.listdir(rundir))

    if tree_name!=None:
        chain = ROOT.TChain(tree_name)
        file_refs.append(chain)
    else:
        chain = ROOT.TChain()

    runs_not_found=list(run_numbers_list[:])

    n_pulses = 0


    for runnumber in run_numbers_list:
        if max_pulses and n_pulses != None:
            if n_pulses > max_pulses:
                break

        for dir in directories:
            dir=rundir+dir+"/"

            assert os.path.isdir(dir)

            abs_path=dir+TTreeBaseName+str(runnumber)+".root"


            if os.path.isfile(abs_path):
                if tree_name==None:
                    _file_ = ROOT.TFile(abs_path)

                    for key in _file_.GetListOfKeys():
                        if key.GetClassName() == "TTree":
                            tree_name = key.GetName()
                            break
                    else:
                        assert False, "could not find tree name from file {0}. \n '_file_.Print() : {1}".format(TTreeBaseName+str(runnumber)+".root",_file_.ls())
                    chain = ROOT.TChain(tree_name)
                    file_refs.append(chain)
                    _file_.Close()


                runs_not_found.remove(runnumber)
                """Calling chain.GetEntries() keeps ROOT from crashing when a non-default value for
                <nentries = TTree::kMaxEntries> is used in TChain::Add. Must be a bug."""
                chain.GetEntries()
                if runnumber in max_events_dict:
                    assert os.path.exists(abs_path)
                    chain.Add(abs_path,max_events_dict[runnumber])
                else:
                    chain.Add(abs_path)

                    if TTreeBaseName == "2ncorr" and n_pulses != None:
                        chain.GetEntry(chain.GetEntries()-1)
                        n_pulses+=chain.EventNumber
                break


        else:
            failed_to_find=map(lambda x:dir+TTreeBaseName+str(x)+".root",runs_not_found)
            failed_to_find="\n".join(failed_to_find)
            assert 0, "\nFailed to locate the following ROOT file{1}: \n {0}".format(failed_to_find,"" if len(runs_not_found)==1 else "s")

    assert isinstance(chain,ROOT.TTree), "Can't find tree in root file \n".format(chain.ls())
    if TTreeBaseName == "2ncorr":
        return chain,n_pulses
    else:
        return chain

class twoNCorrPlotter():
    hist_name_index=1
    def __init__(self,TTree,Min,Max,binwidth):
        self.min=Min
        self.max=Max
        self.nbins=int(float(Max-Min)/binwidth)
        self.binwidth=binwidth

        assert isinstance(TTree,ROOT.TTree,), "argument must be a ROOT TTree instance. "
        self.tree=TTree

    def get_ToF_str(self,angle,offset=0,extra_formatted_cut=""):
        assert angle!=330 and angle != 30, "use get_ToF_bot_str or get_ToF_top_str for segmented detectors. "
        draw_str="0.5(D{0}T+D{0}B)-trig+{1})".format(angle,offset)
        _cut=self.auto_coincidence(draw_str)
        if extra_formatted_cut:
            _cut=cut_AND(_cut,extra_formatted_cut.format(angle,ToF=draw_str))

        return(draw_str,_cut)

    def get_ToF_top_str(self, angle, offset=0):
        draw_str = "D{0}T-trig+{1}".format(angle, offset)
        _cut = self.auto_coincidence(draw_str)
        return (draw_str, _cut)

    def get_ToF_bot_str(self, angle, offset=0):
        draw_str = "D{0}B-trig+{1}".format(angle, offset)
        _cut = self.auto_coincidence(draw_str)
        return (draw_str, _cut)

    def get_ToF_hist(self,angle,offset=0,extra_formatted_cut=""):
        assert angle != 330 and angle != 30, "use get_ToF_bot_hist or get_ToF_top_hist for segmented detectors. "
        _hist=self.get_empty_hist()

        _draw_str,_cut=self.get_ToF_str(angle,offset,extra_formatted_cut)

        self.tree.Project(_hist.GetName(),_draw_str,_cut)

        return _hist

    def get_ToF_bot_hist(self, angle, offset=0):
        _hist = self.get_empty_hist()

        _draw_str, _cut = self.get_ToF_bot_str(angle, offset)

        self.tree.Project(_hist.GetName(), _draw_str, _cut)

        return _hist

    def get_diff_str(self,angle,offset=0,extra_formatted_cut=""):
        draw_str="(D{0}T-D{0}B)+{1}".format(angle,offset)
        _cut=self.auto_coincidence(draw_str)
        if extra_formatted_cut:
            _cut=cut_AND(_cut,extra_formatted_cut.format(angle,ToF="0.5*(D{0}T+D{0}B)-trig+{1}".format(angle,offset)))
        return (draw_str,_cut)

    def get_Diff_hist(self,angle,offset=0,extra_formatted_cut=""):
        _hist=self.get_empty_hist(angle,[-30,30])
        _draw_str,_cut=self.get_diff_str(angle,offset,extra_formatted_cut)
        self.tree.Project(_hist.GetName(),_draw_str,_cut)
        return _hist


    def get_ToF_top_hist(self, angle, offset=0):
        _hist = self.get_empty_hist()

        _draw_str, _cut = self.get_ToF_top_str(angle, offset)

        self.tree.Project(_hist.GetName(), _draw_str, _cut)

        return _hist

    def get_empty_hist(self,angle,range=None):
        if range:
            bin_args=int(float(range[-1]-range[0])/self.binwidth), range[0], range[-1]
        else:
            bin_args=self.nbins, self.min, self.max
        _hist = ROOT.TH1F(str(twoNCorrPlotter.hist_name_index),str(angle),
                          *bin_args)
        twoNCorrPlotter.hist_name_index += 1
        return _hist

    def auto_coincidence(self,draw_str):
        coinc_entries = set(re.findall(r"(D[0-9]+[TB]|trig)", draw_str))
        return cut_allNoneZero(coinc_entries)




"""
Module contents:
    -Inputdecks: Generate an MCNP deck that contain sections of python code with the inputdeck class.
    -Python function to rotate a vector around an arbitrary axes (rotate_matrix())
    -Generate rotation portions of MCNP TRCL card with rotate_matrix
    -Convert Mathmatica expression into python with convert().
    -Fuction ptrac_conv() runs ParseNip on a PTRAC file and organises the directory.
    - Function getplotcosine(theta,phi) will return a string that can be used to orient the MCplot window
      from sperical coordinates.
"""



"""
Class inputdeck(path/to/file,excludelines=,writedir=,numtasks=1)
    YOU MUST COP PAST THIS CLASS INTO YOUR PROGRAM TO RUN IT.
    This class loads the provided MCNP unput deck and evaluates all expressoins enclosed by @'s whithin the variable scope of
    the python session that the class is called form.
    NOTE: avode from module import*, since all variable from this module are now included in the global variable scope.
    <excludelines=> is a list of lines that will be ignored.
    <writedir> This variable specifies a directory to write the MCNP input deck txt file to. Default is the same directory
               of the file, path/to/.
    <excludelines> This is a list of lines to ignor.
    <numtasks> This value will create n identical input decks except the SEED card will be incremented by 2 for each file.
                If numtasks>1, then a bash script named 'cmd.sh' is created in the dir of the input deck. cmd.sh runs all
                of the input decks on separate terminals, automatically.
    Example Usage:
        import numpy as np
        x=6
        y=12
        ***Copy paste class code into your file here ***
        inputdeck("/Users/jeffburggraf/MY_MCNP/Simulations/NeutronPinball2/pinball2.txt")

    Now, every block of code between enclosing @'s will be evaluated as python code within your current scope.
    In this example, a new input deck titled py_pinball2.txt will be created in the ../NeutronPinball2 directory.


    Running a paralell MCNP simulation:
    Set the 'numtasks=' keyword to something greater than 1, and also be sure your deck has this in it: "RAND SEED=@seed@".
    Now, the code will automatically increment the seed and create n input decks. It will also generate a bash file
    named cmd.sh, that when ran, will execute all n input decks.

    The ptrac2ROOT program can concatenate all of the ptrac files via the parseall() function, given that you run it in the same
    directory as all of your ptra* files.
"""




"""Used to generate ptrac names for inputdeck class"""
class alphabet_Increment():
    def __init__(self, n=0, base=26):
        self.N = n
        self.base = base

    def __convert_to_TwentyFive_base__(self, N):
        i = 0
        base_list = []
        if N == 0: return [0]

        while N != 0:
            left_over = N % (self.base ** (i + 1))
            N = N - left_over
            base_list.append(left_over / (self.base ** i))

            i += 1
        base_list.reverse()
        return base_list

    def get_abc(self):
        s = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".lower()
        self.N += 1
        return "".join(map(lambda i: s[i],
                           self.__convert_to_TwentyFive_base__(self.N - 1)))



class inputDeck():
    """
    A tool for using python expressions in MCNP input decks.

    Usage:
        # The write directory defaults as "path/to/".
        # Suppose the input deck contains the following line:
        #   1 so @r*inch@

        inp=inputdeck("path/to/inpdeck",write_dir="same_as_original")

        inch=2.54

        # The code below generate 9 inp files, 3 files with different RND seed's for each of the 3 radii.
        for r in range(3):
            i.write_inp_in_scope(globals(),numthreads=3)
    """

    Inp_number = 0

    Letters=alphabet_Increment(2)

    CommandBlocks=[]

    def __init__(self, path_to_original, write_dir="same_as_original"):
        self.globals = {} # This attr is set when self.write_inp_in_scope() is ivoked.

        assert os.path.isfile(path_to_original), "Input deck doesn't exist."
        self.original_inp_path = path_to_original

        if write_dir == "same_as_original":
            write_dir = os.path.dirname(path_to_original)
        else:
            assert os.path.isdir(
                write_dir), "The 'write_dir' provided is not a directory: write_dir='{0}'".format(
                write_dir)

        if not write_dir[-1] == "/":
            write_dir += "/"

        self.write_dir = write_dir

        self.original_inp_name = os.path.basename(path_to_original)

        self.INP_files_written_by_class=[]

    def CleanDirectory(self):

        for filename in os.listdir(self.write_dir):
            m1 = re.match(r"^py_[0-9]+",filename)
            m2 = re.match(r"autoMultithreadLog[0-9]+\.txt", filename)

            if m1 or m2:
                os.remove(self.write_dir+filename)


        return 1

    @staticmethod
    def __write_cmd_file__(instance, max_osx=1):

        assert isinstance(instance, inputDeck)

        bash_cmds = open(instance.write_dir + "cmd", "w")
        bash_cmds.write("".join(inputDeck.CommandBlocks))
        bash_cmds.close()

        if max_osx:
            os.system("chmod +x {0}cmd".format(instance.write_dir))

        return  1

    def __get_next_new_inp_fileObject__(self,name="default", commandKwrags=None, MACOSX=True,
                                    POLIMI=True):

        next_new_inp_name = "py_{0}{1}".\
            format(inputDeck.Inp_number,self.original_inp_name if name=="default" else name)

        next_new_inp_path = "{0}{1}".format(self.write_dir,
                                            next_new_inp_name)

        self.INP_files_written_by_class.append(next_new_inp_path)

        inputDeck.Inp_number += 1

        letters = inputDeck.Letters.get_abc()

        outpname="out" + letters
        if name!="default":
            outpname+=name

        command = "{5}; {0} i={1} ptrac={2} outp={3} runtpe={4}".format(
            "polimi" if POLIMI else "mcnp6", next_new_inp_name
            , "ptra" + letters, outpname, "runtp" + letters,
            "cd " + self.write_dir)

        if MACOSX:
            BashLine = """osascript -e 'tell app "Terminal"\ndo script "{0} {1} 2>autoMultithreadLog{2}.txt;exit"\nend tell'\n""".format(
                command,
                commandKwrags,inputDeck.Inp_number)
        else:
            BashLine = 'gnome-terminal -x sh -c "!!;{0} {1}"\n'.format(command,
                                                                       commandKwrags)

        inputDeck.CommandBlocks.append(BashLine)
        if len(inputDeck.CommandBlocks) >= 2:
            inputDeck.__write_cmd_file__(self, MACOSX)
        return open(next_new_inp_path, "w")

    def __auto_multiThread__(self,numtasks,globals,newname="default",polimi=True,MAC_OSX=True,Terminal_kwargs=""):

        found_seed_line=0
        temp_inp_path=self.original_inp_path + "temp"

        with open(self.original_inp_path) as file, open(temp_inp_path, "w") as temp:
            blank_line_delimeter_count=0

            if polimi:
                for line in file:

                    if re.match(r"^ *\n *$", line.lower()):
                        blank_line_delimeter_count += 1

                    m = re.match(r"^ *dbcn +.+", line.lower())

                    if m:
                        temp.write("DBCN @seed@\n")
                        found_seed_line = 1
                    else:
                        if blank_line_delimeter_count <= 2:
                            temp.write(line)
                        else:
                            break

                if not found_seed_line:
                    temp.write("DBCN @seed@\n")
            else:

                for line in file:

                    if re.match(r"^ *\n *$", line.lower()):
                        blank_line_delimeter_count += 1

                    m = re.match(r"^ *rand +(?:seed=.+ *|stride=.+)", line.lower())
                    if m:
                        temp.write("RAND SEED=@seed@\n")
                        found_seed_line = 1
                    else:
                        if blank_line_delimeter_count <= 2:
                            temp.write(line)
                        else:
                            break

                if not found_seed_line:
                    temp.write("RAND SEED=@seed@\n")

        temp=open(temp_inp_path)

        for i in range(numtasks):
            globals["seed"] = 2 * i + 1



            new_inp_fileObject = inputDeck.__get_next_new_inp_fileObject__(self,name=newname,
                                                                           commandKwrags=Terminal_kwargs,
                                                                           MACOSX=MAC_OSX,
                                                                           POLIMI=polimi)

            write_inputdeck(temp, new_inp_fileObject, Globals=globals)

        os.remove(temp_inp_path)


        print "Use the commands below to run {} identical simulations simultaneously, all with a different random number seed. \n".format(numtasks)
        print "cd {} \n./cmd\n".format(self.write_dir)

    def write_inp_in_scope(self, Globals, POLIMI=True, MACOSX=True,
                           new_name="default", numthreads=1,
                           MCNPTerminalKwrags=""):
        """
        Create and write an input deck(s). Expressions inside the original file
         that appear between enclosing @'s are evaluated as python code in the
         scope of Globals. If multiple files are created (numthread>1),then  a
         bash file named 'cmd' will automatically run all input decks
         simultaneously in separate shells.

        :param Globals: Global variable scope when evaluating python code.
        :param POLIMI: True for POLIMI, False from MCNP
        :param MACOSX: True for mac, False for Linux
        :param new_name: name of the new input deck. Defauls: "py_n"+"originalName"
        :param numthreads: Automatically make multiple copies of input deck and automatically set random number cards.
        :param MCNPTerminalKwrags: e.g. outp=aaacd .
        :return: None
        """

        self.globals = Globals

        if numthreads==1:
            next_new_inp_file = inputDeck.__get_next_new_inp_fileObject__(self,
                                                                          name=new_name,
                                                                          commandKwrags=MCNPTerminalKwrags,
                                                                          MACOSX=MACOSX,
                                                                          POLIMI=POLIMI)

            with open(self.original_inp_path) as py_input_deck:
                write_inputdeck(py_input_deck, next_new_inp_file, self.globals)
        else :
            print "\nRun the following command to execute all inputfiles:\n\nbash {0}{1}\n".format(self.write_dir,"cmd")
            self.__auto_multiThread__(numthreads,Globals,newname=new_name,polimi=POLIMI,MAC_OSX=MACOSX,
                                      Terminal_kwargs=MCNPTerminalKwrags)

        return  "{0} i={1}  ip".format("polimi" if POLIMI else "mcnp6",self.INP_files_written_by_class[-1])




def write_inputdeck(original_file_object, new_file_object, Globals={}):
    """
    This does the evaluation and writing of input decks passed by the inputdeck class defined above.

    :param path_to_original:  Path to an MCNP input deck in which python code may exist between enclosing @'s.
    :param path_to_new:  Path to a new input deck. The code between enclosing @'s in the new input deck will be
    evaluated as python code within the scope of Globals.
    :param Globals: dictionary of variables and values. To use the scope of your current project, do Globals=globals().
    :return: None
    """
    assert isinstance(new_file_object,file) , "'new_file_object' must be a file instance"
    assert isinstance(original_file_object, file), "'original_file_object' must be a file instance"
    assert isinstance(Globals,dict)

    oldfile=original_file_object # just for brevity.

    oldlines=[line for line in oldfile]

    linenum=1

    newlines=[]

    for line in oldlines: # check for odd number of @'s
        assert len(re.findall(r"@",line))%2==0, "Unmatching @'s on line {0}: {1}".format(linenum,line)


        m=re.finditer(r"@(.*?)@",line)

        line_comment_ReGeX=r"(?i) *c +.+"
        if not re.match(line_comment_ReGeX,line):
            # Commented line are ignored. Otherwise replace all @*@'s in line with evaluated code.
            for match in m:
                pythonCode=match.group(1)
                try:
                    evaluated_result=str(eval(match.group(1),Globals))
                except:
                    evaluated_result = "  ERROR-->({})<--ERROR  ".format(pythonCode)
                    warnings.warn("\nCouldn't evaluate '{0}' in line {1}:\n{2}\n\n{3}\n".format(pythonCode,linenum,line,sys.exc_info()))

                # Use match.group(0) to replace the entire capture with the eval'd code, including the @'s
                line = re.sub(re.escape(match.group(0)), evaluated_result,line)

        newlines.append(line)
        linenum += 1

    new_file_object.write(breaklines3("".join(newlines)))
    new_file_object.seek(0)
    original_file_object.seek(0)


"""Apply to mcnp input deck string. This will fix any cards of length > 80 characters.  """
def breaklines3(string):

    lines=string.split("\n")
    fixedlines=[]
    for line in lines:
        assert isinstance(line,str)
        _splitline=line.split("$")
        line=_splitline[0]

        line=re.sub(r"^(.+[^ ]) +",r"\1 ",line) #remove space at end of line
        line = re.sub(r"^ +$", r"\n", line) # Deal with lines only containing spaces
        if len(_splitline)>=2:
            comment="$"+"$".join(_splitline[1:])
            comment=re.sub(r"(^.*[^ ]) +",r"\1 ",comment)#remove space at end of line
        else:
            comment=""

        if len(line)<=77:
            fixedlines.append("{0}{1}".format(line,comment if comment else ""))
            continue

        if re.match(r"^ {0,4}c",line):
            # for comment lines
            line=re.sub(r"(^ {0,4}c[^ ]*) +",r"\1 ",line)#remove space at end of line
            fixedlines.append(line)
            continue

        s=re.findall(r"(.{1,74})(?: |$)",line)

        if line[0:len(s[0])]!=s[0]:
            missed_portion=line[0:re.search(s[0],line).start()]
            s.insert(0,missed_portion)
            warnings.warn("\nCould not break the following line on a whitespace character: \n\n{}".format(line))

        for block in s:
            if len(block)>=78:
                assert False, "Could not break the following line: "+block
        s=s+[""]
        s=filter(lambda x:x,s)  # Remove any blank strings that snuck in.
        fixedline = "\n     ".join(s)
        fixedline+=comment
        fixedlines.append(fixedline)

    result="\n".join(fixedlines)
    result=re.sub(r"\n( *\n *){1,5}\n", r"\n\n", result,count=2) # remoive tripple lines
    return result



"""
The function convert(<string>) converts mathematica FortranForm[] output into python redable code.
To do this, apply Mathematica's FortranForm[...] function to an expression, then copy the result
as plan text, and pass it to the convert function defined below.
The List() function below must be defined within scope before applying eval()
 to output of the convert() function. """
def List(*stuff):
        return [i for i in stuff]


def convert(string):   #converts the above string into python readable code
    out=re.sub(re.escape("\[Theta]"),"theta",string)
    out=re.sub(r"\n +- +","",out)
    out=re.sub(r" +","",out)
    out=re.sub(r"Cos\(","np.cos(",out)
    out=re.sub(r"Sin\(","np.sin(",out)
    out=re.sub(r"Sqrt\(","np.sqrt(",out)
    out=re.sub(r"E([^a-zA-Z0-9])",r"np.e\1",out)
    return out

"""Below is the Fortran code for a rotation matrix around an arbitrary direction."""
fortrancode="""  List(List((x**4 + x**2*y**2 + y**2*Cos(\[Theta]) + x**2*z**2*Cos(\[Theta]))/(x**2 + y**2),
     -   x*y - x*y*Cos(\[Theta]) - z*Sin(\[Theta]),2*Sin(\[Theta]/2.)*(y*Cos(\[Theta]/2.) + x*z*Sin(\[Theta]/2.))),
     -  List(x*y - x*y*Cos(\[Theta]) + z*Sin(\[Theta]),
     -   (x**2*y**2 + y**4 + x**2*Cos(\[Theta]) + y**2*z**2*Cos(\[Theta]))/(x**2 + y**2),
     -   2*Sin(\[Theta]/2.)*(-(x*Cos(\[Theta]/2.)) + y*z*Sin(\[Theta]/2.))),
     -  List(2*Sin(\[Theta]/2.)*(-(y*Cos(\[Theta]/2.)) + x*z*Sin(\[Theta]/2.)),
     -   2*Sin(\[Theta]/2.)*(x*Cos(\[Theta]/2.) + y*z*Sin(\[Theta]/2.)),z**2 + (x**2 + y**2)*Cos(\[Theta])))"""


"""
TRCL(translationVec=[0,0,0],rotationAxis=[0,0,1],theta=0)
Generate an MCNP TRCL string for use with the cell keyword 'trcl='.
"""
def apply_nested(thing,func):
    assert hasattr(func,"__call__")

    result=[]
    if hasattr(thing,"__iter__"):
        for thing2 in thing:
            result.append(apply_nested(thing2,func))
        return result
    else:
        return func(thing)

def round__(number,ndigits=3):
    _stringQ=False
    ndigits+=1
    if isinstance(number,str):
        number=eval(number)
        _stringQ=True

    rndedNumber=(0.1**ndigits)*(number//(0.1**ndigits))
    rndedNumber=round(rndedNumber,ndigits-1)

    if rndedNumber == 0:
        rndedNumber = 0

    if _stringQ:
        return str(rndedNumber)
    else:
        return rndedNumber


class TRCL():
    """
    TRCL(translationVec=[0,0,0],rotationAxis=[0,0,1],theta=0)
    Generate an MCNP TRCL string for use with the cell keyword 'trcl='.
    """
    def __init__(self,translation=[0,0,0],Axis=[0,0,1],theta=0,round=3):
        self.trans=translation
        self.axis=Axis
        self.theta=theta
        self.matrix=rotate_matrix(self.axis,self.theta,False,rnd=round)
        self.round=round

    def __roundvalue__(self,value):
        if value==0:return 0
        elif value==1:return 1
        else: return round(value,self.round)

    def getCard(self,parentheses=True):
        translationsStr=" ".join([str(i) for i in apply_nested(self.trans,lambda x:round__(x,4))])
        matrix=apply_nested(self.matrix,self.__roundvalue__)
        MCNP_matrix_str= " ".join([" ".join(map(str,a)) for a in matrix])

        if parentheses:
            return "({0} {1})". format(translationsStr,MCNP_matrix_str)
        else:
            return "{0} {1}".format(translationsStr, MCNP_matrix_str)

    def __mul__(self, other):
        assert isinstance(other,TRCL)
        result=TRCL()
        result.matrix=np.dot(self.matrix,other.matrix)
        result.trans=[x1+x2 for x1,x2 in zip(self.trans,other.trans)]
        return result


    def __str__(self):
        return self.getCard()


"""
rotate_matrix(axis,radians,mcnp=,eps=,rnd=)
    This function produces a rotation matrix around an arbitrary axes (by using the Fortran code above).
    <axis> is the axis of rotation aa a simple python list. Doesn't need to be normalized.
    <radians> angle of rotation in radians.
    <mcnp=True> This option has the function return a string that can be used on an MCNP trcl card.
    <eps>=0.00001 This is needed to avoid 0/0 NaN.
    <rnd> Number of digits to round the result to.

applyrotation(matrix,vector)
    Using a matrix returned from rotate_matrix(...), applyrotation(matrix,vector) operates <matrix> on <vector>,
    and returns the result.
"""



def rotate_matrix(axis,radians,mcnp=False,eps=0.0000001,rnd=4):
    x,y,z=axis[0],axis[1],axis[2]

    def List(*stuff):
        return [i for i in stuff]

    norm=np.sqrt(x**2+y**2+z**2)
    stringmatrix=convert(fortrancode)
    matrix=re.sub(r"x","("+str(float(x)/norm)+"+eps"+")",stringmatrix)
    matrix=re.sub(r"y","("+str(float(y)/norm)+"+eps"+")",matrix)
    matrix=re.sub(r"z","("+str(float(z)/norm)+"+eps"+")",matrix)
    matrix=re.sub(r"theta",str(radians),matrix)

    matrix=compile(matrix,"<string>","eval")

    out=eval(matrix)
    rnd=rnd+1
    if mcnp==False:
        Return=[]
        for row in out:
            dummyrow=[]
            for entry in row:

                value = (entry // (0.1 ** rnd)) * (0.1 ** rnd)
                value=round(value,rnd-1)
                dummyrow.append(value)
            Return.append(dummyrow)

    else:
        Return=""
        for row in out:
            dummyrow=""
            for entry in row:
                value = (entry // (0.1 ** rnd)) * (0.1 ** rnd)
                value=round(value,rnd-1)
                dummyrow+=(str(value))
                dummyrow+=" "
            Return+=(dummyrow)
            Return+=" "

    return Return



def applyrotation(matrix,vector,rnd=3):
    result=[]
    for row in matrix:
        result.append(sum(round(a*b,rnd) for a,b in zip(row,vector)))
    return result



def flatten_and_concatenate(*things):
    """
    Construct a single 1-D list containing every non-iterable element nested at
    any level in *args.
    :param things: Can be anything what so ever.
    :return: list of all things.
    """
    result=[]
    if len(things)>1:
        for Obj in things:
            result.extend(flatten_and_concatenate(Obj))
        return result
    else:
        Obj=things[0]
        if hasattr(Obj,"__iter__"):
            for obj in Obj:
                result.extend(flatten_and_concatenate(obj))
            return result
        else:
            return [Obj]


__sample_input_deck__="""c
c
@cell.all_cells()@

c
c
1 so 1.5
2 RPP -150 150 -150 150  -150 150
c

c
SDEF
nps 10
M1000  1001 1  """

"""
# example usage:
import mytools
from mytools import cell,cell_collection
import numpy as np

f=open("__sample_inp_deck__","w")
f.write(mytools.__sample_input_deck__)
f.close()

# Any keyword argument can be used below, e.g. trcl=(...): this would add a trcl entry to the card.
initial_ball=cell(material=1000,geom="-1",rho=-1,imp="imp:n,p=1 ",comment=" Base cell: all others are copied from this with the like_but method. ")

# create a cell_collection for the first set of balls.
collection1=cell_collection()

# add initial_ball to collection1.
collection1.add(initial_ball)

#  The likebut() method creates a duplicate of the cell 'ball',
#  then modifies the duplicate according to the **kwargs added.
#
#  I recommend that when supplying a trcl keyword, use a TRCL instance as the value.
#  However supplying a string will also do.
#
#  The 'add_to_collections' argument is just a convenient way to do the equivalent of:
#     _cell=initial_ball.likebut(trcl=trcl)
#     collection1.add( _cell)
#
for theta in np.arange(0,2*np.pi,2*np.pi/6.0):
    trans_vec=[0,10*np.cos(theta),10*np.sin(theta)]
    trcl=mytools.TRCL(translation=trans_vec,Axis=[1,0,0],theta=theta)

    initial_ball.likebut(add_to_collections=[collection1],trcl=trcl)

# Create another collection named all_balls, initialed with all cells in collection1
# This is needed when defining the geometry of the room in MCNP.
all_balls=cell_collection(collection1)

# To be used for the second iteration in this fractal-style geometry.
collection2=cell_collection()


# On the line where cell.likebut is invoked in the for loop below, notice how 2
# TRCL instances are multiplied. The result is the two transformations are
# performed in series. i.e. the rotation matrices are multiplied and the translation vectors are added.


for theta in np.arange(2*np.pi/6.0/2,2*np.pi,2*np.pi/6.0):
    trans_vec=[0,40*np.cos(theta),40*np.sin(theta)]
    trcl=mytools.TRCL(translation=trans_vec,Axis=[1,0,0],theta=theta)

    for _cell in collection1:
        assert  isinstance(_cell,cell)
        _cell.likebut(add_to_collections=[all_balls,collection2],trcl=_cell.trcl*trcl)

# collection2 is now a hexagonal collection of a hexagonal collections of spheres.
# The loop below creates a hexagon collection of hexagonal collections of hexagonal
# collections of spheres
for theta in np.arange(0,2*np.pi,2*np.pi/6.0):
    trans_vec = [0,66 * np.cos(theta), 65 * np.sin(theta)]
    trcl = mytools.TRCL(translation=trans_vec, Axis=[1, 0, 0], theta=theta)

    for _cell in collection2:
        assert isinstance(_cell, cell)
        _cell.likebut(add_to_collections=[all_balls],
                      trcl=_cell.trcl * trcl)

# The room can be defined usef the cell.intersect() method and the all_balls
# collection. The cell geometry below would read 'the room is the intersection
# of the compliment of all cells in all_balls and the negative sense of surface 2.
room=cell(material=0,geom=cell.intersect(all_balls,"-2"))

# setting imp to None automatically makes the importance of this cell, for all
# particles used thus far, equal to zero.
universe=cell(geom="2",imp=None)

print cell.all_cells()

inp=mytools.inputDeck("/Users/jeffburggraf/PycharmProjects/CellClass/__sample_inp_deck__")
print inp.write_inp_in_scope(globals(),POLIMI=0)
"""

class cell():
    """

    """
    __cellnumber__ = 250
    __cellNumbersUsed__ = set()
    __instances__ = []


    @staticmethod
    def all_cells():
        return "\n".join([c.card for c in cell.__instances__])

    @staticmethod
    def union(*args):
        return "(" + ":".join(cell.__geomParse__(*args)) + ")"

    @staticmethod
    def intersect(*args):
        return "(" + " ".join(cell.__geomParse__(*args)) + ")"

    @staticmethod
    def compliment(*cells_or_collections):
        return "("+ " ".join(["#"+str(cell.cell_number) for cell in flatten_and_concatenate(cells_or_collections)] )+")"

    @staticmethod
    def __geomParse__(*args):
        entries = []
        args=flatten_and_concatenate(args)

        for thing in args:
            if isinstance(thing, cell):
                entries.append("#" + str(thing.cell_number))

            elif isinstance(thing, str):
                entries.append(thing)

            elif isinstance(thing, int):
                entries.append(str(thing))

            elif isinstance(thing, cell_collection):
                for c in thing:
                    entries.append("#" + str(c.cell_number))
            else:
                assert 0, "{0} has invalid type {1} ".format(thing, type(thing))

        return entries

    def __set_kwargs__(self,kwargs,card_entries=[]):
        self.trcl=TRCL()
        for keyword, value in kwargs.iteritems():
            if isinstance(value, TRCL):
                setattr(self, "trcl", value)

                if card_entries:
                    card_entries.append("{0}={1}".format(keyword, value.getCard()))

                continue
            else:
                setattr(self, keyword, value)

            if card_entries:
                card_entries.append("{0}={1}".format(keyword, value))

    def __init__(self,material=0, geom="", rho=-1,imp="imp:n,p=1",comment="", cell_number_override=False, **kwargs):

        self.collections=set([])

        if not cell_number_override:

            while cell.__cellnumber__ in cell.__cellNumbersUsed__:
                cell.__cellnumber__+=1
            self.cell_number = cell.__cellnumber__

        else:
            assert isinstance(cell_number_override,int)
            assert cell_number_override not in cell.__cellNumbersUsed__, "Cannot use cell number {} more than once!".format(cell_number_override)
            self.cell_number  = cell_number_override

        cell.__cellNumbersUsed__.add(self.cell_number)

        while cell.__cellnumber__ in cell.__cellNumbersUsed__:  # set cell.__cellnumber__ to the next unused number.
            cell.__cellnumber__ += 1

        self.geom=geom
        self.material=material
        self.rho=rho
        self.imp=imp
        self.kwargs=kwargs


        card_entries=[]
        card_entries.append(str(self.cell_number))
        card_entries.append(str(material))

        if int(material)!=0:
            card_entries.append(str(rho))

        card_entries.append(self.geom)


        if imp:
            assert isinstance(imp,str)
            card_entries.append(imp)
        else:
            longest=""
            i=0
            for c in cell.__instances__:
                i+=1
                if not c.imp:continue
                if len(c.imp.strip())>len(longest):
                    longest=c.imp.strip()
                if i>20:
                    break
            self.imp=longest[:-1]+"0"
            card_entries.append(self.imp)


        self.__set_kwargs__(kwargs,card_entries)

        if comment:
            card_entries.append("$"+comment)

        self.card=" ".join(card_entries)
        cell.__instances__.append(self)


    def __str__(self):
        return self.card


    def likebut(self,add_to_collections=[],**kwargs):
        assert kwargs, "**kwargs cannot be empty. "
        s="{0} like {1} but ".format(cell.__cellnumber__,self.cell_number)
        s+=" ".join(["{0}={1}".format(keyword,value) for keyword,value in kwargs.iteritems()])

        if "mat" in kwargs:
            new_mat=kwargs["mat"]
            kwargs.pop("mat")
        else:
            new_mat=self.material

        if "rho" in kwargs:
            new_rho=kwargs["rho"]
            kwargs.pop("rho")
        else:
            new_rho=self.rho

        newcell=cell(material=new_mat,geom="",rho=new_rho,imp=self.imp,**kwargs)
        newcell.card=s

        if add_to_collections:
            assert isinstance(add_to_collections,list) , "add_to_collections arg must be list of collections to add to. "
            for collection in add_to_collections:
                assert isinstance(collection,cell_collection)
                collection.add(newcell)

        return newcell



class cell_collection():
    def __init__(self,*cells_or_collections):

        self.cells=flatten_and_concatenate(cells_or_collections)

        assert all([isinstance(c,cell) for c in self.cells ]), "all args passed to cell_collection must be a cell instance"

        for _cell in self.cells:
            _cell.collections.add(self)

        self.tags={}



    def add(self,cells_or_collections,tag="",):
        assert isinstance(tag,str),"Second argument to add ('tag') must be a str. "

        cells=flatten_and_concatenate(cells_or_collections)
        assert all([isinstance(c, cell) for c in cells]), "All elements passed with the first argument of cell_collection.add() must be a cell instance"

        self.cells.extend(cells)

        for _cell in cells:
            _cell.collections.add(self)

        if tag:
            if tag in self.tags:
                self.tags[tag]=self.tags[tag].update(cells)
            else:
                self.tags[tag]=set(cells)

    def get_by_tags(self,exclude_tags=[],include_tags=[]):
        """
        :param exclude_tags:
        :param include_tags:
        :return: cell_collection instance of cells w/ tags which pass the include and exclude lists
        """
        _set = set()

        for tag in include_tags:
            if tag in self.tags:
                _set=_set.update(self.tags[tag])
        for tag in exclude_tags:
            if tag in self.tags:
                _set-=self.tags[tag]

        return cell_collection(_set)

    def cell_numbers_str_list(self):
        return [str(c.cell_number) for c in self.cells]


    def __iter__(self):
        return self.cells[:].__iter__()

    def likebut(self,add_to_collection=None,**kwargs):
        new_cells=[c.likebut(add_to_collections=0,**kwargs) for c in list(self)]

        if isinstance(add_to_collection,cell_collection):
            add_to_collection.add(new_cells)
            return add_to_collection
        else:
            assert add_to_collection==None, "add_to_collection must be either None or a cell_collection instnace. "
            return cell_collection(new_cells)

    def __str__(self):
        return "\n".join([str(c) for c in self])

    def __getitem__(self, item):
        return self.cells[item]





class ptracFile():
    """
    class ptracFile(pathtoptracfile):
        This class has a method for converting a ptrac file with ParsNip.

        METHODS:
            ptrac_conv:
                This method converts ptrac file with ParsNip.
                Example:
                file=ptracFile("/Users/jeffburggraf/MY_MCNP/Simulations/CrossTalk2/ptrac")
                test.ptrac_conv(pathtoptracfile="/Users/jeffburggraf/MY_MCNP/ParsNIP",verbose=1,nLines="3")
                ...> Now you must compy past the commands that print to the screen.
                ...> <pathtoptracfile> is the absolute path to the folder containg the parsnip src folder.
                ...> <nLines> is the number of input lines in the ptrac file. usuallyt three.
    """
    def __init__(self,pathtoptracfile):
        assert  os.path.isfile(pathtoptracfile),"ERROR: ptrac file: "+str(pathtoptracfile)+" does not exist!"
        self.ptracfilepath=pathtoptracfile
        self.outpdir=os.path.dirname(pathtoptracfile)


    def ptrac_conv(self,parseniplocation="/Users/jeffburggraf/MY_MCNP/ParsNIP",gccCompile=True,nLines="3"):
        if parseniplocation[-1]!="/":
            path=parseniplocation+"/"  #path = Parsnip home folder
        else:
            path=parseniplocation

        assert os.path.isdir(path+"src"), "No 'src' directory in the parseniplocation you provided. "
        assert all(map(lambda x: os.path.isfile(x),[path+"src/"+"main.c",path+"src/"+"lookup.c",path+"src/"+"main_h.h"])), "You 'src' folder must contain the folloing files: 'main.c', 'lookup.c', and main_h.h'"
        self.nipdir=self.outpdir+"/ParsNip/"
        if not os.path.isdir(self.nipdir):
            answer=raw_input("Must create a sub-directory titled 'ParseNip' whithin "+self.outpdir+". (y/n)")
            if answer=="y" or answer=="Y" or answer=="yes":
                os.mkdir(self.nipdir)
            else:
                assert False


        dummycomandfile=open(self.nipdir+"ParsNipCommands.sh","w")
        dummycomandfile.write("#!/bin/bash\n")

        if gccCompile:
            print "\nPasting the commands below into a terminal will convert the ptrac files with Parsnip:\n"
            print "cd "+path
            print "gcc src/*.c -o parsnip"
            print "cd "+self.nipdir
            print "chmod u+x ParsNipCommands.sh"
            print "./ParsNipCommands.sh\n"
        else:
            print "\nPasting the commands below into a terminal will convert the ptrac to ROOT:\n"
            print "cd "+self.nipdir
            print 'root -l'
            print 'const char* str="'+self.nipdir+'Parsenip_inp"'
            print ".L "+path+"src/main.c"
            print "run(str)"

        inp_template="#\n"+path+"\n#\nptrac_in\n#\nptrac_out\n#\n1\n#\n#\n"+nLines+"\n#\n#\n0"
        dummy_parse_input=[i for i in inp_template.splitlines()]
        dummy_parse_input[3]=self.ptracfilepath
        dummy_parse_input[5]="ptrac_out"
        parse_input=open(self.nipdir+"Parsenip_inp","w")
        for line in dummy_parse_input:
            parse_input.write(line+"\n")
        parse_input.close()
        dummycommands= "cd "+self.nipdir+"\n"
        dummycommands+=path+"parsnip "+"Parsenip_inp"+"\n"
        for line in dummycommands.splitlines():
            dummycomandfile.write(line+"\n")
        dummycomandfile.close()

"""
getplotcosine(theta,phi)
    This function will return a string that can be used to orient the MCplot window
    from spherical coordinates with the BASIS command.
"""
phihat="""        List(-((Sin(phi)*Sin(eps + theta))/
     -     Sqrt(Abs(Cos(phi)*Sin(eps + theta))**2 + Abs(Sin(phi)*Sin(eps + theta))**2)),
     -  (Cos(phi)*Sin(eps + theta))/
     -   Sqrt(Abs(Cos(phi)*Sin(eps + theta))**2 + Abs(Sin(phi)*Sin(eps + theta))**2),0)"""

def getplotcosine(theta,phi,eps=0.001):
    def List(*stuff):
        return [i for i in stuff]
    def Abs(num):
        return np.abs(num)
    horizontal=eval(convert(phihat))
    verticle=[np.cos(phi)*np.cos(theta),np.cos(theta)*np.sin(phi),np.sin(theta)]
    return "BASIS  "+" ".join([str(round(value,3)) for value in horizontal+verticle])


"""
Used to Draw a list of drawable root objects. Finds an optimal way to split 45 or less object among canvii.
Example usage:

    histos=[]
for i in range(13):
    histos.append(ROOT.TH1F(str(i),str(i),100,0,10))
    histos[-1].FillRandom("gaus",10000);

Drawer=mytools.drawAllRootHistos(histos)
Drawer.Draw("hist")
"""

class drawAllRootHistos(object):
    def __init__(self,histlist,width=1000):
        N=len(histlist)
        self.histos=histlist
        self.w=width;self.h=int(self.w*9/16.0)
        if "c1" in globals():
            self.basename="C"
        else:
            self.basename="c"
        self.c={i:None for i in range(1,6)}
        PadcdNums={i:None for i in range(1,5+1)}
        if N==1:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            PadcdNums[1]=[1]
        elif N==2:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[1].Divide(2);PadcdNums[1]=[1,2]
        elif 3<=N<=4:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[1].Divide(2,2);PadcdNums[1]=[1,2,3,4]
        elif 5<=N<=8:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[1].Divide(2,2);PadcdNums[1]=[1,2,3,4]
            self.c[2].Divide(2,2);PadcdNums[2]=[1,2,3,4]
        elif N==9:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
        elif N==10:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
            PadcdNums[2]=[1]
        elif 11<=N<=13:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
            self.c[2].Divide(2,2);PadcdNums[2]=[1,2,3,4]
        elif 14<=N<=18:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
            self.c[2].Divide(3,3);PadcdNums[2]=range(1,9+1)
        elif 19<=N<=27:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[3]=ROOT.TCanvas(self.basename+"3", self.basename+"3", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
            self.c[2].Divide(3,3);PadcdNums[2]=range(1,9+1)
            self.c[3].Divide(3,3);PadcdNums[3]=range(1,9+1)
        elif 28<=N<=9*4:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[3]=ROOT.TCanvas(self.basename+"3", self.basename+"3", self.w, self.h)
            self.c[4]=ROOT.TCanvas(self.basename+"4", self.basename+"4", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
            self.c[2].Divide(3,3);PadcdNums[2]=range(1,9+1)
            self.c[3].Divide(3,3);PadcdNums[3]=range(1,9+1)
            self.c[4].Divide(3,3);PadcdNums[4]=range(1,9+1)
        elif (9*4+1)<=N<=9*5:
            self.c[1]=ROOT.TCanvas(self.basename+"1", self.basename+"1", self.w, self.h)
            self.c[2]=ROOT.TCanvas(self.basename+"2", self.basename+"2", self.w, self.h)
            self.c[3]=ROOT.TCanvas(self.basename+"3", self.basename+"3", self.w, self.h)
            self.c[4]=ROOT.TCanvas(self.basename+"4", self.basename+"4", self.w, self.h)
            self.c[5]=ROOT.TCanvas(self.basename+"5", self.basename+"5", self.w, self.h)
            self.c[1].Divide(3,3);PadcdNums[1]=range(1,9+1)
            self.c[2].Divide(3,3);PadcdNums[2]=range(1,9+1)
            self.c[3].Divide(3,3);PadcdNums[3]=range(1,9+1)
            self.c[4].Divide(3,3);PadcdNums[4]=range(1,9+1)
            self.c[5].Divide(3,3);PadcdNums[5]=range(1,9+1)
        else:
            assert False, "Too many histograms to Draw!"
        self.PadcdNums=PadcdNums


    def Draw(self,Drawoptions=""):
        histnum=0;
        print self.PadcdNums
        for canvasindex,cdindexlist in self.PadcdNums.iteritems():
            if cdindexlist:
                for cdindex in cdindexlist:
                    if histnum==len(self.histos):
                        return 1
                    self.c[canvasindex].cd(cdindex)
                    self.histos[histnum].Draw(Drawoptions)
                    histnum+=1



"""
Iterator for histogram bins:
for bin in histIter(hist):
    bin.binnum  # current bin number
    bin.x  # bin center
    bin.y # bin content
    bin.error # bin error
    bin.xl  # lower edge of bin
    bin.width #

"""

class histIter():
    class bin():
        def __init__(self,num,x,xleft,y,error,width,numbins,hist):
            self.binnum = num
            self.x = x
            self.xleft = xleft
            self.y = y
            self.error = error
            self.width = width
            self.numbins = numbins
            self.hist=hist

    def __init__(self,hist):
        self.binnum=0
        self.width = float(hist.GetBinWidth(1))
        self.numbins=hist.GetNbinsX()
        self.hist=hist

        self.overflow=self.hist.GetBinContent(self.numbins+1)
        self.underflow = self.hist.GetBinContent(0)

    def __iter__(self):
        return self

    def next(self):
        self.binnum += 1
        if self.binnum == 1 + self.numbins:
            raise StopIteration
        x=self.hist.GetXaxis().GetBinCenter(self.binnum)
        xl = x-self.width/2.0

        y=self.hist.GetBinContent(self.binnum)
        error=self.hist.GetBinError(self.binnum)
        return histIter.bin(self.binnum,x,xl,y,error,self.width,self.numbins,self.hist)

"""
HistToList gives the raw data from a ROOT histogram. x values are bin centers.
"""
class HistToList(object):
    def __init__(self,hist):
        if isinstance(hist,ROOT.TH1):

            assert isinstance(hist,ROOT.TH1F)
            length=hist.GetNbinsX()
            self.hist=hist
            self.n=length
            self.x=np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1,length+1)],np.float32)
            self.y=np.array([hist.GetBinContent(i)for i in range(1,length+1)],np.float32)
            self.xyTuple=[(hist.GetXaxis().GetBinCenter(i),hist.GetBinContent(i)) for i in range(1,length+1)]
            self.hist=hist
            self.erry=np.array([hist.GetBinError(i)for i in range(1,length+1)],np.float32)
            self.width = hist.GetBinWidth(1)
            self.binLeftEdges = self.x-self.width/2.0

            non_zero_indicies = np.where(self.y>0)
            self.x_NonZero = self.x[non_zero_indicies]
            self.y_NonZero = self.y[non_zero_indicies]
            self.erry_NonZero = self.erry[non_zero_indicies]

            self.binnum=1

    def MaximumBinError(self,x):
        if self.hist.GetBinLowEdge(1) <= x < (self.hist.GetBinLowEdge(self.n)+self.width):
            return self.hist.GetBinError(1+np.argmax(self.x>=x))
        else:
            return 0

    def __iter__(self):
        return histIter(self.hist)

    def __len__(self):
        return self.n

class HistToList2D:
    def __init__(self,hist,x_axis="x",y_axis="y"):
        x_axis = eval("hist.Get{}axis()".format(x_axis.upper()))
        y_axis = eval("hist.Get{}axis()".format(y_axis.upper()))
        self.xyValues = np.zeros((x_axis.GetNbins(), y_axis.GetNbins()), dtype=np.float32)
        self.xyErrors = np.zeros((x_axis.GetNbins(), y_axis.GetNbins()), dtype=np.float32)
        self.xCenters = np.zeros(x_axis.GetNbins(), dtype=np.float32)
        self.yCenters = np.zeros(y_axis.GetNbins(), dtype=np.float32)

        for ix in range(0, x_axis.GetNbins()):
            for iy in range(0, y_axis.GetNbins() ):
                self.xyValues[ix][iy] = hist.GetBinContent(ix+1, iy+1)
                self.xyErrors[ix][iy] = hist.GetBinError(ix+1, iy+1)

                self.xCenters[ix] = x_axis.GetBinCenter(ix+1)
                self.yCenters[iy] = y_axis.GetBinCenter(iy+1)

        self.xyCenters = zip(self.xCenters, self.yCenters)

class HistToList3D:
    def __init__(self, hist, x_axis="x", y_axis="y", z_axis="z"):

        if  isinstance(x_axis,str):
            x_axis = eval("hist.Get{}axis()".format(x_axis.upper()))
        else:
            assert isinstance(x_axis,ROOT.TAxis)

        if  isinstance(y_axis,str):
            y_axis = eval("hist.Get{}axis()".format(y_axis.upper()))
        else:
            assert isinstance(y_axis, ROOT.TAxis)

        if  isinstance(z_axis,str):
            z_axis = eval("hist.Get{}axis()".format(z_axis.upper()))
        else:
            assert isinstance(z_axis, ROOT.TAxis)

        self.xyzValues = np.zeros((x_axis.GetNbins(), y_axis.GetNbins(),z_axis.GetNbins()), dtype=np.float32)
        self.xyzErrors = np.zeros_like(self.xyzValues, dtype=np.float32)
        self.xyzVolumes = np.zeros_like(self.xyzValues)

        self.xCenters = np.zeros(x_axis.GetNbins(), dtype=np.float32)
        self.xLeftEdges = np.zeros(x_axis.GetNbins(), dtype=np.float32)
        self.xRightEdges = np.zeros(x_axis.GetNbins(), dtype=np.float32)

        self.yCenters = np.zeros(y_axis.GetNbins(), dtype=np.float32)
        self.yLeftEdges = np.zeros(y_axis.GetNbins(), dtype=np.float32)
        self.yRightEdges = np.zeros(y_axis.GetNbins(), dtype=np.float32)

        self.zCenters = np.zeros(z_axis.GetNbins(), dtype=np.float32)
        self.zLeftEdges = np.zeros(z_axis.GetNbins(), dtype=np.float32)
        self.zRightEdges = np.zeros(z_axis.GetNbins(), dtype=np.float32)

        self.xBinWidths = np.zeros_like(self.xCenters)
        self.yBinWidths = np.zeros_like(self.yCenters)
        self.zBinWidths = np.zeros_like(self.zCenters)


        for ix in range(0, x_axis.GetNbins()):
            for iy in range(0, y_axis.GetNbins()):
                for iz in range(0, z_axis.GetNbins()):
                    self.xyzValues[ix][iy][iz] = hist.GetBinContent(ix + 1, iy + 1, iz + 1)
                    self.xyzErrors[ix][iy][iz] = hist.GetBinError(ix + 1, iy + 1, iz + 1)

                    self.xCenters[ix] = x_axis.GetBinCenter(ix + 1)
                    self.xLeftEdges[ix] = x_axis.GetBinLowEdge(ix+1)
                    self.xRightEdges[ix] = x_axis.GetBinUpEdge(ix+1)

                    self.yCenters[iy] = y_axis.GetBinCenter(iy + 1)
                    self.yLeftEdges[iy] =y_axis.GetBinLowEdge(iy + 1)
                    self.yRightEdges[iy] =y_axis.GetBinUpEdge(iy + 1)

                    self.zCenters[iz] = z_axis.GetBinCenter(iz + 1)
                    self.zLeftEdges[iz] = z_axis.GetBinLowEdge(iz + 1)
                    self.zRightEdges[iz] = z_axis.GetBinUpEdge(iz + 1)

                    self.xBinWidths[ix] =  x_axis.GetBinWidth(ix+1)
                    self.yBinWidths[iy] =  y_axis.GetBinWidth(iy+1)
                    self.zBinWidths[iz] =  z_axis.GetBinWidth(iz+1)

                    self.xyzVolumes[ix][iy][iz] = self.xBinWidths[ix]*self.yBinWidths[iy]*self.zBinWidths[iz]


        self.xyzCenters = zip(self.xCenters, self.yCenters, self.zCenters)

class Legend(object):
    def __init__(self,histlist=None, fontsize=0.05,TCanvas=None,xmin=0.78, ymin=0.7, xmax=1, ymax=0.95,title=""):

        self.legEntriesArgs=[]
        self.NEntries = 0
        self.fontsize=fontsize;self.TCanvas=None;self.xmin=xmin;self.ymin=ymin;self.xmax=xmax;self.ymax=ymax
        self.xmin=xmin;self.xmax=xmax;self.ymin=ymin;self.ymax=ymax



        if histlist:
            self.leg = ROOT.TLegend(xmin, ymin, xmax, ymax)
            if title:
                self.setTitle(title)
            self.leg.SetTextSize(fontsize)
            self.TCanvas = TCanvas

            for HistogramHistTitle in histlist:
                self.addentry(HistogramHistTitle[0], HistogramHistTitle[1], "f")

            self.leg.Draw()

    def Draw(self,Drawoption="",fontsize=False,TCanvas=None, xmin=0.76, ymin=0.8, xmax=0.89, ymax=0.89,title=""):
        self.leg = ROOT.TLegend(xmin, ymin, xmax, ymax)
        if self.NEntries==0:return 0
        self.leg.SetLineColorAlpha(ROOT.kBlack,0)
        if title:
            self.setTitle(title)
        if fontsize:
            self.leg.SetTextSize(fontsize)
        else:

            self.leg.SetTextSize(self.fontsize)

        for args in self.legEntriesArgs:
            self.leg.AddEntry(*args)

        self.leg.Draw(Drawoption)
        if self.TCanvas:
            if isinstance(self.TCanvas, ROOT.TCanvas):
                self.TCanvas.cd()

    def addentry(self,hist,legEntry,AddEntryoption=""):
        self.NEntries+=1
        if isinstance(hist,ROOT.TH1):AddEntryoption+="f"
        if isinstance(hist, ROOT.TGraph): AddEntryoption += "l"
        if isinstance(hist, ROOT.TH2F): AddEntryoption += "P"
        self.legEntriesArgs.append( (hist, str(legEntry), AddEntryoption) )

    def setTitle(self,title):
        self.leg.SetHeader(str(title))



def ListToHist(x,y,title="auto",binCenters_NotLeftEdge=True):
    assert len(x) >= 2, "Must have more than a single bin"

    bin_widths = [x[i]-x[i-1] for i in range(1,len(x))]

    assert all([dx>0 for dx in bin_widths]) , "x list must be in ascending order. "

    assert len(x)==len(y) ," x and y lists must be of same length. "

    avg_bin_width=np.average(bin_widths)
    assert all([abs(avg_bin_width-width)/float(avg_bin_width)<0.01 for width in bin_widths])  ,  "All bins must be of nearly (+- 1%) equal length, i.e. x[i+1]-x[i]==c+-0.01 for all i. "


    title=str(title)
    if title=="auto":
        if hasattr(ListToHist,"auto_title_numbers"):
            ListToHist.auto_title_numbers+=1
        else:
            ListToHist.auto_title_numbers=0
        title="hist{}".format( ListToHist.auto_title_numbers)

    nbins=len(x)

    if binCenters_NotLeftEdge:
        Xmin=x[0]-avg_bin_width/2.0
        Xmax = x[-1] + avg_bin_width / 2.0
    else:
        Xmin=x[0]
        Xmax=x[-1]+avg_bin_width

    _hist=ROOT.TH1F(title,title,nbins,Xmin,Xmax)

    [_hist.SetBinContent(i+1,y[i]) for i in range(nbins)]

    return _hist


def ListToHistOld(Xmin,Xmax,y,title="hist0"):
    nbins=len(y)

    if title=="hist0":
        # increment the auto-title
        if hasattr(ListToHist,"titlesUsed"):
            ListToHist.titlesUsed+=1
        else:
            ListToHist.titlesUsed=0
        title="hist{}".format(ListToHist.titlesUsed)

    Hist=ROOT.TH1F(title,title,nbins,Xmin,Xmax)

    for i in range(nbins):
        Hist.SetBinContent(i + 1,y[i]);

    return Hist

from numpy.polynomial import legendre as leg
import scipy.integrate as integrate


class LegendreSeries():
    """
    LegendreSeries(TH1F hist, int n, bool LeastSquares) ; Computes the Legendre series expansion to order n.
    self.TF1: root function of the series result.

    If LeastSquares==True, then the histogram is stripped of any zero bins, and a least squares fit to a Legendre
    polynomial of order n is performed. Bin errors are assumed to be Poissonian for the fit.
    Use of the LeastSquares method is preferred when there are many empty bins, or i.e. when the histogram is "gappy".
    """
    __num_calls__ = 0  # for ROOT TF1 names

    def __init__(self,hist, MaxCoeff_OR_CoeffList = 3,LeastSquares=False):

        assert isinstance(hist,ROOT.TH1)

        hist_func = lambda x: hist.Interpolate(x)

        if hasattr(MaxCoeff_OR_CoeffList,"__iter__"):
            coeffs_included = np.array(MaxCoeff_OR_CoeffList)
        elif isinstance(MaxCoeff_OR_CoeffList,int):
            coeffs_included = range(MaxCoeff_OR_CoeffList+1)
        else:
            assert 0, "MaxCoeff_OR_CoeffList must be a list of ints, or a single integer equal to the max order to be used.  "

        if LeastSquares:
            hist_list = HistToList(hist)
            x,y, erry = hist_list.x_NonZero, hist_list.y_NonZero, hist_list.erry_NonZero
            self.coeff = leg.legfit(x,y, coeffs_included, w= 1.0/erry)
        else:
            self.coeff = []
            for c in coeffs_included:
                _b = legendre(c)
                result = integrate.quad(lambda x: _b(x)*hist_func(x), -1, 1)
                self.coeff.append((c+0.5)*result[0])


        def _func(x, *pars):
            return leg.legval(x[0], self.coeff)

        self.TF1 = ROOT.TF1("LefExpansion"+str(LegendreSeries.__num_calls__),_func, -1,1,0)
        self.TF1.SetLineColor(hist.GetLineColor())

        LegendreSeries.__num_calls__+=1
"""
Fit a histogram with legendre poly's up to 10nth degree.
!!NOTE: THIS FUNCTION WORKS IN THE COSINE OF THETA!!
Histogram must be filled with direction cosine only.
The 'self.fit' method returns the TF1 object.

Example: q=LegeneraHistFit(hist)
         q.fit(4,evenonly=True,FitStatistics=False)

self.fit method arguments:
        degree: The max order of  legendre polynomial to use
        evenonly: Option to use only even-order poly's
        functionname: Function name at the ROOT end.
        FitStatistics: Whether of not to draw the fit stats box.

Attributes:
    self.hist: Histogram
    self.parameters: Fit parameters in ascending order.

"""
class LegeneraHistFit():
    def __init__(self,hist):
        assert isinstance(hist,(ROOT.TH1,ROOT.TGraph)), "Error: Argument must be TH1 or TGraph object"
        ROOT.gSystem.Load("libMathMore");
        self.hist=hist

    def fit(self,degree,evenonly=False,functionname="LegFit",FitStatistics=True):
        def generateLegFuncString(values):
            return "+".join(["[{0}]*ROOT::Math::legendre({0},x)".format(i) for i in values])
        if evenonly:
            coefs=[2*k for k in range(degree)]
        else:
            coefs=[k for k in range(degree)]
        assert len(coefs)<=10, "Error: Too many parameters for TF1."
        self.fitfunc=ROOT.TF1(functionname,generateLegFuncString(coefs),-1,1)
        self.fitfunc.SetParameters(*tuple([1]*len(coefs)))
        self.fitfunc.SetParNames(*tuple(["P_{"+"{}".format(c)+"}(x)" for c in coefs]))
        self.hist.Fit(self.fitfunc)
        self.fitfunc.Draw("same")
        self.parameters=np.frombuffer(self.fitfunc.GetParameters(),count=len(coefs))
        if FitStatistics:
            ROOT.gStyle.SetOptFit(1011)
        return self.fitfunc



"""
Three dimensional histogram class, ThreeDimHist.
Example:

    h=ROOT.TH3F("","",  100,0,10,  10,0,10,  10,0,10)
    h=ThreeDimHist(10,0,10, 10,0,10,  10,0,10,Titlefunc=lambda x: x+" [cm]")
    for i in range(10000):
        z=(10*np.random.uniform()**2)
        h.Fill((np.random.randn()+5,np.random.randn()+5,z))
    hist=h.GetZSliceHist(3)
    hist.Draw("colz")
"""
class ThreeDimHist():
    def __init__(self,nbinsx,xlow,xhigh,nbinsy,ylow,yhigh,nbinsz,zlow,zhigh,Titlefunc=lambda i:str(i)):
        assert zlow<zhigh and ylow<yhigh and xlow<xhigh, "[z,y,z]high must be greater than its corresponding low!"
        self.zbinwidth=float(zhigh-zlow)/float(nbinsz)
        self.xbinwidth=float(xhigh-xlow)/float(nbinsx)
        self.ybinwidth=float(yhigh-ylow)/float(nbinsy)
        self.binvolume=self.xbinwidth*self.ybinwidth*self.zbinwidth
        self.nbinz=nbinsz;self.zlow=zlow;self.zhigh=zhigh;

        def gettitel(binnum):
            if binnum==-1:
                return "Underflow"
            elif binnum>=nbinsz:
                return "Overflow"
            else:
                return str(zlow+self.zbinwidth*binnum)
        self.histograms={i:ROOT.TH2F(Titlefunc(gettitel(i)),Titlefunc(gettitel(i)),nbinsx,xlow,xhigh,nbinsy,ylow,yhigh) for i in range(-1,nbinsz+1)}


    def Fill(self,(x,y,z),weight=1):
        bin=int(float(z-self.zlow)//self.zbinwidth)
        if 0<=bin<=(self.nbinz-1):
            self.histograms[bin].Fill(x,y,weight)
        elif bin>=self.nbinz:
            self.histograms[self.nbinz].Fill(x,y,weight)
        else:
            self.histograms[-1].Fill(x,y,weight)

    def GetZSliceHist(self,z):
        bin=int(float(z-self.zlow)//self.zbinwidth)
        if 0<=bin<=(self.nbinz-1):
            return self.histograms[bin]
        elif bin>=self.nbinz:
            return self.histograms[self.nbinz]
        else:
            return self.histograms[-1]

    def GetBinHistIter(self):
        # doesn't include over/underflow
        return [(str(self.zlow+self.zbinwidth*int(key)),value) for key,value in list(self.histograms.iteritems())[0:-2]]




"""
Gaussian broden from list. Example:
    h=gaussianbroaden(range(100),map(lambda x: 0.98*x,range(100)),sigma=2)
    h.original.Draw(" hist")
    h.original.Set
    h.result.Draw("hist same")
    # Can also supple a histogram as argument
    hh=gaussianbroaden(hist=h.original, sigma=1)
    c2=ROOT.TCanvas()
    hh.result.Draw("hist")
"""

class gaussianbroaden(object):
    def __init__(self,x=None,y=None,sigma=1,hist=False):
        if x and y:
            n=len(x)
            assert n==len(y), "x and y must be of equal length. "
            assert all([x[i+1]>x[i] for i in range(n-1)]) ,"x list must be monotonically increasing. "
            integral=sum(y)
            min=x[0];max=x[-1]
            binwidth=float(max-min)/n
        elif hist:
            assert isinstance(hist,ROOT.TH1F),  " 'hist' must be a TH1F object. "
            n=hist.GetNbinsX()
            binwidth=hist.GetBinWidth(1)
            max=hist.GetMaximumBin()+binwidth
            min=hist.GetMinimumBin()-binwidth
            integral=hist.Integral()
            HH=HistToList(hist)
            x=HH.x;y=HH.y

        normconstant=1./(np.sqrt(2.*np.pi)*float(sigma))
        centerdeviation=int(4*sigma//binwidth)
        inithist=ROOT.TH1F("dummy","dummy",n,min-centerdeviation*binwidth,max+centerdeviation*binwidth)
        [inithist.Fill(float(xx),float(yy)) for xx,yy in zip(x,y)]
        bincenters=[inithist.GetXaxis().GetBinCenter(i) for i in range(1,n+1)]
        binheigths=[inithist.GetBinContent(i)for i in range(1,n+1)]
        def gaussian(x,mean): return normconstant*np.e**(-float((mean-x)**2)/(2.*sigma**2.))
        def GetConvolutedBinHeight(binindex,binvalue):
                bincenter=bincenters[binindex]
                if centerdeviation==0 or binvalue==0:
                    return [(binindex,binvalue)]
                return map(lambda x: (x,binvalue*gaussian(x,bincenter)),bincenters[binindex-centerdeviation:binindex+centerdeviation])
        resulthist=ROOT.TH1F("dummy2","dumm2",n,min-centerdeviation*binwidth,max+centerdeviation*binwidth)
        for index,center in enumerate(binheigths):
            [resulthist.Fill(*_) for _ in GetConvolutedBinHeight(index,center)]
        self.original=inithist
        self.original.SetLineColor(ROOT.kOrange)
        self.result=resulthist

"""
Convolute a root TH1.
"""

def ConvolveHist(hist,sigma,overwritehist=True,newhistTitle="_histnew"):
    assert isinstance(hist,ROOT.TH1)
    if not overwritehist:
        hist=hist.Clone(newhistTitle)
    integral=float(hist.Integral())
    histtolist=HistToList(hist)
    binwidth=histtolist.width
    sigma=sigma/float(binwidth)
    kernal=signal.gaussian(int(6*sigma), std=sigma)
    result=ndimage.filters.convolve1d(histtolist.y, kernal)
    for index,bincontent in enumerate(result):
        index+=1
        hist.SetBinContent(index,bincontent)
    renorm=integral/float(hist.Integral())
    hist.Scale(renorm)
    return hist


def vectorAngle(v1,v2,radians=True,eps=0.00001):
    dot=sum([x1*x2 for x1,x2 in zip(v1,v2)])
    mag1mag2=np.sqrt(sum([x*x for x in v1]))*np.sqrt(sum([x*x for x in v2]))

    cos=dot/(mag1mag2)
    if cos>1.-eps or cos<-1.+eps:
        return 0
    else:
        return np.arccos(cos) if radians else (180/np.pi)*np.arccos(cos)



def subsets(iterable, length,__level__=0,__q__=None):
    """
    :param iterable: set
    :param length: length of subsets to be generated
    :param __level__: don't use
    :param __q__: don't use
    :return: A list containing all possible subsets of length <length> contained in <iterable>.
    """
    self=subsets
    if __q__==None:

        self.result=[]
        self.q=np.zeros(length,dtype=np.int)
        self.N=len(iterable)

    if __level__==length:
        self.result.append(itemgetter(*self.q.copy())(iterable))
    else:
        start=self.q[__level__-1]+1 if __level__>0 else 0
        for _x in range(start,self.N):
            self.q[__level__]=_x
            subsets(iterable,length,__level__+1, 1)

        return  np.array(self.result)

"""
Homemade vector class akin to the Vpython vector class. Vectors support all stardard vector operations.
Examples:
    vec1=vector((1,1,1))
    vec2=vector((0,1,0)).unit(len=2)  # Resize the vector. Deflt is len=1

    vec1,vec2
    >>>(1, 1, 1) (0.0, 2.0, 0.0)

    vec2%vec1                       #Cross product
    >>>(2.0, 0.0, -2.0)

    vec2/2                          #Divide by constant
    (0.,1.,0.)

    vec2.scale(10)                  #Multiply by constant
    >>>(0.0, 20.0, 0.0)

    vec2+vec2                       #Add (subtract also works).
    >>>(1.0, 3.0, 1.0)

    abs(vec1-vec2)                  #Magnitude
    >>>1.73

    vec2**3                       #|vec2|**3
    >>>8.0

    vec1.sphericalform(degrees=False,ndigits=4)       #Convert to spherical form [r,theta,phi]
    >>>[1.7321, 0.9553, 0.7854]

    vec2.cylindricalform(degrees=False,ndigits=4)    #Convert to cyndircal form [r,theta,z]
    >>>[2.0, 1.5708, 0.0]

    vec1.mathematicaform()
    >>>{1, 1, 1}

    vec3=vec1-vec2.unit().scale((vec1*vec2.unit()))  #Gram-Schmidt orthogonalization
    >>>(1.0, 0.0, 1.0)
    vec3*vec2
    >>>0.0

    -vec2                                            #Negative
    >>>(-0.0, -2.0, -0.0)

    vec2==--vec2                                     # Equality
    >>>True

    vec1.x,vec1.y,vec1.z                            # Attributiues
    >>>1,1,1

    vector((1,np.pi/4,np.pi/4),coordinates="spherical")  #Specify in spherical coordinates. May also use cylindrical.
    >>>(0.500, 0.499, 0.707)

    vec1.angle(vec2,degree=False)                   #Vector angle
    >>>0.955
"""

class vector():
    def __init__(self,sequence,coordinates="cartesian",degrees=False):
        if coordinates=="cartesian":
            self.array=np.array(sequence)
            self.x=self.array[0]
            self.y=self.array[1]
            if len(sequence)==3:
                self.z=self.array[2]
            else:
                self.z=1
            self.current=0

        elif coordinates=="spherical":
            if degrees:
                k=np.pi/180
            else:
                k=1
            r,th,phi=tuple(sequence)
            th=k*th;phi=k*phi
            self.x=r*np.sin(th)*np.cos(phi)
            self.y=r*np.sin(th)*np.sin(phi)
            self.z=r*np.cos(th)
            self.array=np.array((self.x,self.y,self.z))
        elif coordinates=="cylindrical":
            if degrees:
                k=np.pi/180
            else:
                k=1
            r,phi,z=tuple(sequence)
            phi=k*phi
            self.x=r*np.cos(phi);self.y=r*np.sin(phi);self.z=z
            self.array=np.array((self.x,self.y,self.z))
        else:
            assert False, "ERROR:'"+str(coordinates)+"'"+' is not a valid coordinate system!\n\t\t\t    Options are "cartesian","spherical",or "cylindrical" specified ' \
                          'by (x,y,z) (r,theta,phi) and (r,phi,z), respectively. '

    def __repr__(self):
        return str(tuple(self.array))

    def __str__(self):
        return self.__repr__()


    def __add__(self, other):
        return vector([x+y for x,y in izip_longest(self.array,other)])

    def __sub__(self, other):
        return vector([x-y for x,y in izip_longest(self.array,other)])

    def __abs__(self):
        return np.sqrt(sum([x**2 for x in self.array]))

    def __mul__(self, other):
        return sum([x*y for x,y in izip_longest(self.array,other)])

    def __iter__(self):
        return self

    def __mod__(self, other):
        u1,u2,u3=tuple(self.array);v1,v2,v3=tuple(other)
        return vector((u2*v3-u3*v2,u3*v1-u1*v3,u1*v2-u2*v1))

    def __getitem__(self, item):
        return self.array[item]

    def unit(self,len=1):
        mag=abs(self)
        return vector([len*x/mag for x in self.array])

    def __div__(self, other):
        other=float(other)
        return vector([x/other for x in self.array])

    def __pow__(self, power, modulo=None):
        return abs(self)**power

    def scale(self,k):
        return vector([x*k for x in self.array])

    def mag(self):
        return abs(self)

    def next(self):
        if self.current >2:
            self.current=0
            raise StopIteration
        else:
            self.current += 1
            return self.array[self.current-1]
    def __eq__(self, other):
        return np.allclose([x for x in self],[x for x in other])
    def __neg__(self):
        return vector([-x for x in self.array])

    def angle(self,other,eps=0.00001,degree=True):
        dot=sum([x1*x2 for x1,x2 in izip_longest(self,other)])
        mag1mag2=abs(self)*np.sqrt(sum([x*x for x in other]))
        cos=dot/(mag1mag2)
        if cos>1.-eps or cos<-1.+eps:
            return 0
        elif degree==True:
            return (180/np.pi)*np.arccos(cos)
        else:
            return np.arccos(cos)

    def sphericalform(self,degrees=False,ndigits=4):
        # (r, theta, phi)
        if degrees:
            k=180/np.pi
        else:
            k=1.
        x,y,z=tuple(self.array)
        return map(lambda x:round(x,ndigits),[abs(self),k*np.arctan(np.sqrt(x**2+y**2)/z),k*np.angle(np.complex(real=x,imag=y))])

    def cylindricalform(self,degrees=False,ndigits=4):
        # (r, theta, phi)
        if degrees:
            k=180/np.pi
        else:
            k=1.
        x,y,z=tuple(self.array)
        return  map(lambda x:round(x,ndigits),[np.sqrt(x**2+y**2),k*np.angle(np.complex(real=x,imag=y)),z])

    def mathematicaform(self):
        return re.sub("\)","}",re.sub("\(","{",str(self)))


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.
    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.
    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.
    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])
    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out



"""Iterate through particles of the same nps number.
Use on TTrees generated by Jeff's mppost.c or parse.c (polimiparse.c) program

    Iterator that returns iterators which loop through each particle history in the simulation.

    Example:
    for nps in mytools.treeiter(tree):
       for event in nps:
            pass
            # everything in this block loop has the same particle nps.
            # using break is allowed.
        """


class treeiter():
    class npsiter():
        def __init__(self, treeiter):
            self.tree = treeiter.tree
            self.ientry = treeiter.ientry

            self.tree.GetEntry(self.ientry)
            _nps = self.tree.nps
            while _nps == self.tree.nps:
                treeiter.ientry += 1
                self.tree.GetEntry(treeiter.ientry)
                if treeiter.ientry >= treeiter.max:
                    treeiter.ientry = treeiter.max
                    break
            # At this point, treeiter.ientry is either at the beginning of the next particle history
            # or at the end of the entire tree loop.

            self.max = treeiter.ientry

        def __iter__(self):
            return self

        def __back_an_event__(self):
            self.ientry-=1
            self.tree.GetEntry(self.ientry)

        def next(self):
            if self.ientry >= self.max:
                raise StopIteration

            self.tree.GetEntry(self.ientry)
            self.ientry += 1
            return self.tree

    def __init__(self,tree,maxEntries=None):
        """
        Iterator that returns iterators which loop through each particle history in the simulation.
        \n
        Example:
        for nps in mytools.treeiter(tree):
           for event in nps:
                pass
                # everything in this block loop has the same particle nps.

        """
        assert isinstance(tree,ROOT.TTree)
        self.tree=tree
        self.entries=tree.GetEntries()
        self.ientry=0
        self.max=maxEntries if maxEntries else self.entries
        if self.max>self.entries: self.max=self.entries

    def __iter__(self):
        return self

    def next(self):

        if self.ientry>=self.max:
            raise  StopIteration
        else:
            return treeiter.npsiter(self)

class POLIMI_tree_iter():
    """
    Triple iterator

    Example:

    for nps in POLIMI_tree_iter(tree):
        for src_particle_history in nps:
            for event in src_particle_history:
                # event is a reference to a TTree instance, so all it's branches are attributes of the variable event.
                # the first event here is always a POLIMI src particle birth, whether it's a bnk event or a src event (only the first src particle is a src event in POLIMI).


    Iterator model:
        Each level defines a method called stopcheck(self) which, using the global self.tree and self.ientry[] variables,
        raises StopIteration if such conditions are met for that level. In the beginning of the next method of each iterator, the stopchecker methods
        are call for ALL outer iterators and the current iterator. All next event incriments (GetEntry and ientry+=1) occur in
        the next method of the innermost iterator.

    """
    def __init__(self,tree,Max_entries=None):
        assert isinstance(tree,ROOT.TTree)

        self.MMMaxentry=Max_entries-1 if (Max_entries and Max_entries<=tree.GetEntries()) else tree.GetEntries()-1
        self.tree=tree

        self.ientry=[0]
        self.tree.GetEntry(0)

    def stopcheck(self,Raise=True):
        if self.ientry[0] > self.MMMaxentry:
            if Raise:
                raise StopIteration
            else:
                return 1
        else:
            return 0

    def __iter__(self):
        return self

    def next(self):
        self.stopcheck()
        return POLIMI_tree_iter.npsiter(self)

    class npsiter():
        def __init__(self, treeiter):
            assert isinstance(treeiter, POLIMI_tree_iter)
            self.mother = treeiter
            self.tree = self.mother.tree
            self.ientry = self.mother.ientry
            self.initial_nps = self.tree.nps

            # assure that we are at the beginning of an nps, if not go to the next one.
            if self.ientry[0]!=0:
                self.tree.GetEntry(self.ientry[0]-1) #  look at pervious event
                if self.tree.nps==self.initial_nps:
                    while self.tree.nps==self.initial_nps:
                        self.ientry[0]+=1
                        self.tree.GetEntry(self.ientry[0])
                        # print self.ientry[0],self.initial_nps
                        if self.mother.stopcheck(Raise=False):
                            break
                #     we now must be at the beginning of the next nps
                    self.initial_nps = self.tree.nps
                else:
                    self.tree.GetEntry(self.ientry[0])

            _has_run=0
            while self.tree.event == 1:
                _has_run=1
                self.ientry[0] += 1
                self.tree.GetEntry(self.ientry[0])
                if self.stopcheck(Raise=0) or self.mother.stopcheck(Raise=0):
                    break
            if _has_run:
                self.ientry[0] -=1

                # print self.ientry[0]
            self.tree.GetEntry(self.ientry[0])


        def stopcheck(self,Raise=1):
            if self.tree.nps != self.initial_nps:
                if Raise:
                    raise StopIteration
                else:
                    return 1
            else:
                return 0

        def __iter__(self):
            return self

        def next(self):
            self.mother.stopcheck()
            self.stopcheck()
            return POLIMI_tree_iter.npsiter.srciter(self)

        class srciter():
            def __init__(self, npsiter):
                assert isinstance(npsiter,POLIMI_tree_iter.npsiter)
                self.mother=npsiter
                self.grandmother=self.mother.mother
                self.ientry=self.mother.ientry
                self.tree=self.mother.tree

                self.next_calls=0

                while self.tree.bnknum!=25 and self.tree.bnknum!=26 and self.tree.event!=1:
                    self.next_calls+=1
                    self.ientry[0]+=1
                    self.tree.GetEntry(self.ientry[0])
                    if self.mother.stopcheck(Raise=0) or self.grandmother.stopcheck(Raise=0) or self.stopcheck(Raise=0):
                        break


            def stopcheck(self,Raise=1):
                if self.next_calls and (
                        self.tree.bnknum == 25 or self.tree.bnknum == 26):
                    if Raise:raise StopIteration
                    else: return 1
                else:
                    return 0


            def __iter__(self):
                return self

            def next(self):
                self.tree.GetEntry(self.ientry[0])

                self.grandmother.stopcheck()
                self.mother.stopcheck()



                self.stopcheck()

                self.next_calls += 1

                self.ientry[0]+=1

                return self.tree


channel_dict = {(1, 1): "D30T",   (1, 2): "D30B",
                (1, 3): "D54T",   (1, 4): "D54B",
                (1, 5): "D78T",   (1, 6): "D78B",
                (1, 7): "D102T",  (1, 8): "D102B",
                (1, 9): "D126T", (1, 10): "D126B",
                (1, 11):"D150T", (1, 12): "D150B",
                (1, 13):"D210T", (1, 14): "D210B",
                (1, 15): "D234T", (1, 16): "D234B",
                (1, 17): "D258T", (1, 18): "D258B",
                (1, 19): "D282T", (1, 20): "D282B",
                (1, 21): "D306T", (1, 22): "D306B",
                (1, 23): "D330T", (1, 24): "D330B",
                }

class TDC1190_DXXTB_to_DrawStr_str():
    def __init__(self,runnumber):
        tree=NCorrTree(runnumber,TTreeBaseName="r")

        _branch_names_=[]
        found=0

        for _branch in tree.GetListOfBranches():
            for _leaf in _branch.GetListOfLeaves():
                _branch_names_.append(_leaf.GetName())
                if "1190" in _branch_names_[-1]:
                    found=1
                    break
            if found:
                break
        else:
            assert False, "Cannot find TDC1190 leaf. The names of all leaves in TTree are: {0}".format(_branch_names_)


        self.TDC1190_name=_branch_names_[-1]

    def get_draw_str(self,angle,top_Or_bot="top",hit_num=1):
        angle=str(angle)

        angles=set([re.match(r"D([0-9]+)",_v_).group(1)  for _v_ in channel_dict.values()])
        assert angle in angles, "Invalid detector angle, {}".format(angle)
        top_Or_bot=top_Or_bot.lower()
        assert top_Or_bot in ["t","b","top","bot"], " Argument top_ORbot must be one of the following: 't' 'b','top','bot'"
        topbot="T"if top_Or_bot in ["t","top"] else "B"

        ch=channel_dict.keys()[channel_dict.values().index("D{0}{1}".format(angle,topbot))][-1]
        return "{0}[{1}][{2}]".format(self.TDC1190_name,hit_num,ch)

def get_array_index(hit,rootVar):
    """TDC1190[100][129]/I<->TDC1190[hit][rootvar]/I"""
    return hit*129+rootVar



def conv(path_to_rootFile,start_entry=0,stop_entry=None):
    assert os.path.isfile(path_to_rootFile)
    runnumber = int(re.search(r"r([0-9]+)\.root$", path_to_rootFile).groups(0)[0])

    file=ROOT.TFile(path_to_rootFile)
    oldtree = file.Get("R1DC")
    if start_entry:
        oldtree.GetEntry(start_entry)

    newpath = os.path.dirname(path_to_rootFile) + "/newr{0}.root".format(runnumber)

    if start_entry:
        newfile = ROOT.TFile(newpath, "update")
        newtree = newfile.Get("tree")
        newtree.GetEntry(newtree.GetEntries())
        assert oldtree.GetEntries()>newtree.GetEntries()
    else:
        newfile=ROOT.TFile(newpath,"recreate")
        newtree=ROOT.TTree("tree","tree")

    # Dictionary mapping R1DC1190[i][j] to new leaf name, e.g. R1DC1190[1][1]->D30T
    dict = channel_dict = {(1, 1): "D30T", (1, 2): "D30B", (1, 3): "D90T", (1, 4): "D90B",
        (1, 5): "D150T", (1, 6): "D150B", (1, 7): "D210T",
        (1, 8): "D210B", (1, 9): "D120T", (1, 10): "D120B",
        (1, 11): "D270T", (1, 12): "D270B", (1, 13): "D330T",
        (1, 14): "D330B", (1, 16): "trig"}

    addresses=[array.array("f",[0]) for key,value in sorted(dict.iteritems())]
    names=[value for key,value in sorted(dict.iteritems())]
    indicies=[get_array_index(key[0],key[1]) for key,value in sorted(dict.iteritems())]

    for name,address in zip(names,addresses):
        if start_entry:
            newtree.SetBranchAddress(name, address)
        else:
            newtree.Branch(name,address,name+"/F")


    address_indicie=[(a,i) for a,i in zip(addresses,indicies)]

    for i in xrange(start_entry,stop_entry if stop_entry else oldtree.GetEntries()-1):

        for address,index in address_indicie:
            address[0]=oldtree.TDC1190[index]*0.1

        newtree.Fill()

    assert (oldtree.GetEntries()==newtree.GetEntries() or stop_entry!=None), "Entry mismatch between trees: {0}  VS  {1}".format(oldtree.GetEntries(),newtree.GetEntries())
    newtree.Write()



class newTree():
    """
    Creates a new tree with the same branch structure and bindings as a tree.
    Calling the self.write_current_nps() writes all events of the particle
    history at the tree's current read entry. \n

    REMEMBER TO CALL "self.Write" AFTER RUNNING LOOP.

    To apply a cut on events written where the cut is the same through out the entire
    analysis, do the following:
        t=newtree(tree,un_changing_cut="zaid==6000")
        t.write_current_nps_const_cut()

    To apply a cut that may change with each nps, do the following:
        t=newtree(tree)
        t.write_current_nps_changing_cut(changingcut)
    """

    def __init__(self, MotherTTree, DaughterTTreeName="auto", un_changing_cut=None):
        self.oldtree = MotherTTree
        self.max=self.oldtree.GetEntries()
        self.isDone=0 # to prevent double Writing on accident.
        _readentry = MotherTTree.GetReadEntry()

        if DaughterTTreeName == "auto":
            basename = os.path.basename(self.oldtree.GetDirectory().GetName())
            self.filename = "new_" + basename
        else:
            self.filename = DaughterTTreeName

        self.newfile = ROOT.TFile(self.filename, "recreate")
        self.newfile.cd()
        self.newtree = self.oldtree.CloneTree(0)

        self.previousy_written_NPS=-1

        MotherTTree.GetEntry(_readentry)  # Place the tree marker back where it was, just in case
        self.oldtree.GetDirectory().cd() # move back into file

        self.priotNPSWritten=0

        if un_changing_cut:
            self.Tformula = ROOT.TTreeFormula("nps_cut", un_changing_cut, self.oldtree)
            self.Tformula.GetNdata()
        else:
            self.Tformula = 1

    def write_current_nps(self, __cut__=lambda: 1):
        nps = self.oldtree.nps
        if self.priotNPSWritten==nps:
            return 0



        if self.previousy_written_NPS==nps:
            return 0
        else:
            self.previousy_written_NPS = nps

        ientry_initial = self.oldtree.GetReadEntry()
        ientry = ientry_initial


        if ientry > 1:
            while self.oldtree.nps == nps:
                ientry -= 1
                self.oldtree.GetEntry(ientry)
                if ientry<=0: # ientry==0 is the first event of three. ientry's <0 are meaningless.
                    ientry=0
                    break
            else:
                ientry += 1  # Go back to start of particle history

        # Write until the position at the time this function was called
        for ientry in xrange(ientry, ientry_initial + 1):
            self.oldtree.GetEntry(ientry)
            if __cut__():
                self.newtree.Fill()

        ientry = ientry_initial + 1
        if self.oldtree.GetEntry(ientry)==0: # The end of the tree has been reached if ...==0
            return 1

        while self.oldtree.nps == nps:
            if ientry >= self.max:
                break
            if __cut__():
                self.newtree.Fill()

            ientry += 1
            self.oldtree.GetEntry(ientry)

        self.priotNPSWritten=nps
        return 1

    def write_current_nps_changing_cut(self, cutString):
        TTreeFormula = ROOT.TTreeFormula("nps_cut", cutString, self.oldtree)
        self.write_current_nps(cut=lambda: TTreeFormula.EvalInstance())

    def write_current_nps_const_cut(self):
        self.write_current_nps(cut=lambda: self.Tformula.EvalInstance())

    def write_current_event(self):
        self.newtree.Fill()

    def Write(self):
        if self.isDone==1:
            warnings.warn("Tree already written! Cannot invoke newtree.Write() more than once. Function call ignored. ")
            return 0
        self.isDone=1
        self.newfile.cd()
        self.newtree.Write()
        return 1


"""
Equiv. to Mathematica's RandomChoice[].
"""
def weighted_values(values, probabilities, size=1):
    bins = np.add.accumulate(probabilities)
    bins=bins/bins[-1]
    return np.array(values)[np.digitize(np.random.random_sample(size), bins)]

#  Same as above. For when the distribution needs to be sampled from many times.
class weighted_values_class():
    def __init__(self,values, probabilities):
        bins = np.add.accumulate(probabilities)
        self.bins=np.array(bins)/np.array(bins)[-1]
        self.values=np.array(values)
    def getRandom(self,size=1):
        return self.values[np.digitize(np.random.random_sample(size), self.bins)]





#  Replace every key in pattern_replaceDict with the corresponding value. Regex expressions with capturing groups are allowed.
#  This function only replaces the first occurrence from the dictionary. See recursiveReplaceAll() to replace all occurences.
def findfirstAndReturnSplit(pattern_replaceDict,string):
    assert isinstance(pattern_replaceDict,dict), " First argument must be a dictionary. e.g.  {pattern1:replacement1,..}. Regex are allowed."
    _matches=[re.search(pattern,string) for pattern in pattern_replaceDict.keys()]
    dictpositions=[]
    matches=[]
    for i,match in enumerate(_matches):
        if match:
            matches.append(match)
            dictpositions.append(i)
    if not matches:
        return False

    bestMatch=matches[0]
    index=0
    bestMatchIndex=0
    for match in matches[1:]:
        index+=1
        if match.start()>bestMatch.start():continue
        elif match.start()<bestMatch.start():
            bestMatch=match
            bestMatchIndex=index
            continue
        elif match.start()==bestMatch.start():
            if match.end()>bestMatch.end():
                bestMatch=match
                bestMatchIndex=index
            elif match.end()==bestMatch.end():
                print "Warning: double match with equal length at same location: 1: {}  2: {} .  Choose option 1".format(bestMatch.group(),match.group())
            else: continue
    bestMatchIndex=dictpositions[bestMatchIndex]
    leftstop=bestMatch.start()
    rightstop=bestMatch.end()
    stringreplacement=re.sub(pattern_replaceDict.keys()[bestMatchIndex],pattern_replaceDict.values()[bestMatchIndex],string[leftstop:rightstop])
    leftstring=string[:leftstop]+stringreplacement
    rightstring=string[rightstop:]
    return [leftstring,rightstring]

# The first argument is a dictionary with keys as regex expressions that will be replaced with the corresponding values. Capture groups are allowed.
def recursiveReplaceAll(pattern_replaceDict,string):
    result=""
    remaining=string
    while 1:
        _dummy=findfirstAndReturnSplit(pattern_replaceDict,remaining)
        if _dummy:
            remaining=_dummy[1]
            result+=_dummy[0]
        else:
            result+=remaining
            break
    return "".join(result)


"""
Use mathematica's FortranForm[] function on any expression, which may include all of the commonly used math functions, to return a string produce a string that can
 then be parsed into either numpy or root with the function below.
 This function can also round decimals.

"""
def fortranFormToPython(string,roundNdigits=4,language="python"):
    language=language.lower()
    assert isinstance(string,str), "Argument must be string"
    listbounder = ("[", "]") if language.lower() == "python" else ("{", "}")
    string=re.sub(r"\n {2,}- {2,}","",string)
    string = re.sub(r"\n", "", string)
    string = re.sub(r" +", "", string)
    string=string.replace("\n","")
    string=re.sub(r"(\\\[((?i)[a-z]+)\])",lambda m:m.group(2).lower(),string)  # Replace Mathematica symbols with lover case strings
    string=re.sub("Sec\(","1/Cos(",string)
    # pi and E sub
    if language!="python":
        string=re.sub("E\*\*","TMath::E()**",string)
        string=re.sub("Pi", "TMath::Pi()", string)
    else:
        string = re.sub("E\*\*", "np.e**", string)
        string = re.sub("Pi", "np.pi", string)

    specialFunctionKeys=['Cos', 'Sin', 'Tan', 'ArcCos', 'ArcSin', 'ArcTan', 'Cosh', 'Sinh', 'Tanh', 'ArcCosh', 'ArcSinh', 'ArcTanh', 'Log',"Abs","Sqrt"]

    if language=="python":
        specialFunctionValues=["np."+i.lower() for i in specialFunctionKeys]
    else:
        specialFunctionValues=[]
        for key in specialFunctionKeys:
            if key[0:3]=="Arc":
                key="A"+key[3:]
            if key[-1]=="h":
                key=key[0:-1]+"H"
            specialFunctionValues.append("TMath::"+key)

    specialFunctionDict={pattern:repl for pattern,repl in zip(specialFunctionKeys,specialFunctionValues)}
    round_pattern_entry={r"([0-9]+\.[0-9]{1,"+str(roundNdigits)+"})[0-9]+":r"\1"}
    specialFunctionDict.update(round_pattern_entry)

    def ParseFirstList(string):
        match=re.search("List\(",string)
        if not match:
            return False
        rp=1;lp=0
        begin=match.end()
        end=None
        right=listbounder[1]
        left=listbounder[0]
        for pos,char in enumerate(string[begin:]):
            if char=="(":rp+=1
            if char==")": lp+=1
            if rp==lp:
                end=begin+pos
                break
        if end==None:
            print "Warning: unbalanced parentheses in: "+string[begin-5:]
            return False
        else:
            return string[:begin-5]+left+string[begin:end]+right+string[end+1:]
    while True:
        newstring=ParseFirstList(string)
        if newstring==False:
            break
        else:
            string=newstring
    return recursiveReplaceAll(specialFunctionDict,string)


"""
Passing a ROOT TH1 to this function will return a text block containing an SI and SP card pair, which allows MCNP to sample
from the same statistics as the ROOT histogram. The same can be done with a callable python object, provided that the min,max, and nbins
is given.

If histInterpolation=True (default), then a smooth probability distribution is constructed from the histogram.
Otherwise, values are sampled uniformly from a given bin.
"""
def mcnp_probability_dist(Th1_or_callable,SINumber,min=None,max=None,nbins=None,nDigits=4,histInterpolation=True):
    assert SINumber>0, "SInumber must be an integer greater than zero. "
    SINumber=int(SINumber)
    option= "A" if histInterpolation else "H"
    if isinstance(Th1_or_callable,ROOT.TH1):
        hist=Th1_or_callable
        histotolist=HistToList(hist)
        xbins=[x-histotolist.width/2.0 for x in histotolist.x]
        if min:
            _xbins = [x for x in xbins if x >= min]
            _truncated_x_min_position=len(xbins)-len(_xbins)
        else:
            _truncated_x_min_position=0
        if max:
            _xbins = [x for x in xbins if x <= max]
            _truncated_x_max_position=-(len(xbins)-len(_xbins))+1
            if _truncated_x_max_position>=0:
                _truncated_x_max_position=len(xbins)
        else:
            _truncated_x_max_position=len(xbins)

        xbins = xbins[_truncated_x_min_position:_truncated_x_max_position]

        xbins = xbins+[xbins[-1]+histotolist.width]

        ybins=list(histotolist.y[_truncated_x_min_position:_truncated_x_max_position])+[0.0]
        assert nbins==None, "Don't use nbins arg for TH1, just rebin the histogram instead. "

    elif hasattr(Th1_or_callable,"__call__"):
        from scipy import integrate
        assert min!=None and max!=None and nbins!=None, "All of the following keyword arguments must be  supplied: min=, max=, nfuncbin= "
        binwidth=float(max-min)/float(nbins)
        xbins=[min+n*binwidth for n in range(nbins+1)]
        ybins=[0.0]
        for x1,x2 in zip(xbins[:-1],xbins[1:]):
            Sum=integrate.quad(Th1_or_callable,x1,x2)[0]
            ybins.append(Sum)


    else:
        assert 0, "Th1_or_callable must be a ROOT TH1 or __call__ able object. "

    assert len(xbins)==len(ybins), "Error: bin lengths are unequal. There are problems, check code..."
    ybins=np.array(ybins,dtype=np.float64)/float(sum(ybins))
    SI=" ".join([str(round(num,nDigits)) for num in xbins])
    SP=" ".join([str(round(num,nDigits)) for num in ybins])
    result="SI{0} {1} ".format(SINumber,option)+SI+"\n"+"SP{0}   ".format(SINumber)+SP
    return result





def removeExcess(string):
    assert isinstance(string, str)
    allgroups = []
    groups = []

    counter = 0
    for char in string:
        if char == "(":
            counter += 1;
            groups.append("")
            continue
        elif char == ")":
            counter -= 1;
            continue

        for i in range(0, counter):
            groups[i] += char
        if counter == 0 and groups:
            [allgroups.append(g) for g in groups]
            counter = 0
            groups = []

def removeExcess1(string):
    m= re.match(r"([^\(\)]*)\((.+)\)([^\(\)]*)",string)
    # print m.groups()
removeExcess1("21212((a+b)(c)**(x+e^^(14-6))+14)-4")

"\((?:[^()]|(?0))*\)"



class mathString(str):
    def __init__(self, String, remove_paret=1):
        self.__needParentheses=1
        self.string = String

    def __add__(self, other):
        return mathString("{0}+{1}".format(self, other))

    def __radd__(self, other):
        return mathString("{0}+{1}".format(other, self))

    def __sub__(self, other):
        return mathString("{0}-{1}".format(self, other))

    def __rsub__(self, other):
        return mathString("{0}-{1}".format(other, self))

    def __div__(self, other):
        other = mathString(other)
        if self.__needParentheses and other.__needParentheses:
            Format = "({0})/({1})"
        elif self.__needParentheses:
            Format = "({0})/{1}"
        elif other.__needParentheses:
            Format = "{0}/({1})"
        else:
            Format = "{0}/{1}"

        return mathString(Format.format(self, other))

    def __rdiv__(self, other):
        self, other = mathString(other), self
        if self.__needParentheses and other.__needParentheses:
            Format = "({0})/({1})"
        elif self.__needParentheses:
            Format = "({0})/{1}"
        elif other.__needParentheses:
            Format = "{0}/({1})"
        else:
            Format = "{0}/{1}"

        return mathString(Format.format(self, other))

    def __and__(self, other):
        return mathString("({0}) && ({1}) ".format(self, other))

    def __rand__(self, other):
        return mathString("({0}) && ({1}) ".format(other, self))

    def __or__(self, other):
        return mathString("({0}) || ({1}) ".format(self, other))

    def __ror__(self, other):
        return mathString("({0}) || ({1}) ".format(other, self))

    def __gt__(self, other):
        return "{0}>{1}".format(self,other)

    def __lt__(self, other):
        return "{0}<{1}".format(self,other)

    def __pow__(self, power, modulo=None):
        return mathString("({0})**({1})".format(self,power))

    def __mul__(self, other):
        return mathString("({0})*({1})".format(self,other))

    def __rmul__(self, other):
        return self.__mul__(other)

ROOTColors=[4,2,3,1,7,6,41,8,11,5,12,30]


def set_palette(name="palette", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array.array('d', stops)
    r = array.array('d', red)
    g = array.array('d', green)
    b = array.array('d', blue)

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)


def get_array_index(hit,rootVar):
    """TDC1190[100][129]/I<->TDC1190[hit][rootvar]/I   ... Used in R1DC1990()"""
    return hit*129+rootVar



def R1DC1990_conv(path_to_rootFile,start_entry=0,stop_entry=None,channel_map_dict=None,updating_file=False):
    """
    Restructures a daq2 r####.root TTree. Able to map selected R1DC1990[i][j]
    channels to a more user friendly leaf name. Resulting TTree MUCH faster and smaller.
    :param path_to_rootFile: Absolute pth to r####.root file.
    :param start_entry: Update
    :param stop_entry:
    :param channel_map_dict: Dictionary that maps R1DC1190[i][j]'s to new leaf names, e.g. R1DC1190[1][1]->D30T in the example below (See variable "dict" below. )
    :return:
    """
    runnumber = int(
        re.search(r"r([0-9]+)\.root$", path_to_rootFile).groups(0)[0])

    # If root file not foud in given directory, check the download directory and move it to correct directory (symptom of laziness)
    if not os.path.isfile(path_to_rootFile):
        baseRootFileName=os.path.basename(path_to_rootFile)
        downloadsPath="/Users/jeffburggraf/Downloads/"
        if os.path.isfile(downloadsPath+baseRootFileName):

            copiedTreePath=os.path.dirname(path_to_rootFile)+"/"+baseRootFileName
            print downloadsPath+baseRootFileName,copiedTreePath

            copyfile(downloadsPath+baseRootFileName,copiedTreePath)

            assert os.path.isfile(copiedTreePath), "Coipying root file from download directory failed"

            R1DC1990_conv(copiedTreePath,start_entry,stop_entry,channel_map_dict)
            return None

        else:
            assert 0, "file {} cannot be found. ".format(path_to_rootFile)


    file=ROOT.TFile(path_to_rootFile)
    oldtree = file.Get("R1DC")


    newpath = os.path.dirname(path_to_rootFile) + "/newr{0}.root".format(runnumber)

    if updating_file:
        assert os.path.isfile(newpath), "Cannot update {}, since it hasn't been created. "
        newfile = ROOT.TFile(newpath, "update")
        newtree = newfile.Get("tree")
        if oldtree.GetEntries() == newtree.GetEntries():
            return 1
        newtree.GetEntry(newtree.GetEntries())
        assert oldtree.GetEntries()>newtree.GetEntries()
        start_entry=newtree.GetEntries()
        print "Updating run number {0} from {1} entries to {2} entries...".format(runnumber,newtree.GetEntries(),oldtree.GetEntries())
    else:
        newfile=ROOT.TFile(newpath,"recreate")
        newtree=ROOT.TTree("tree","tree")


    if channel_map_dict==None:
        # Dictionary mapping R1DC1190[i][j] to new leaf name, e.g. R1DC1190[1][1]->D30T in the exaple below
        dict =  {(1, 1): "D30T", (1, 2): "D30B", (1, 3): "D90T", (1, 4): "D90B",
        (1, 5): "D150T", (1, 6): "D150B", (1, 7): "D210T",
        (1, 8): "D210B", (1, 9): "D120T", (1, 10): "D120B",
        (1, 11): "D270T", (1, 12): "D270B", (1, 13): "D330T",
        (1, 14): "D330B", (1, 16): "trig"}
    else:
        dict=channel_map_dict

    addresses=[array.array("f",[0]) for key,value in sorted(dict.iteritems())]
    names=[value for key,value in sorted(dict.iteritems())]
    indicies=[get_array_index(key[0],key[1]) for key,value in sorted(dict.iteritems())]

    for name,address in zip(names,addresses):
        if updating_file:
            newtree.SetBranchAddress(name, address)
        else:
            newtree.Branch(name,address,name+"/F")


    address_indicie=[(a,i) for a,i in zip(addresses,indicies)]

    for i in xrange(start_entry,stop_entry if stop_entry else oldtree.GetEntries()):
        oldtree.GetEntry(i)
        for address,index in address_indicie:
            address[0]=oldtree.TDC1190[index]*0.1

        newtree.Fill()


    if (oldtree.GetEntries()!=newtree.GetEntries() or stop_entry!=None):
        warnings.warn( "Entry mismatch between trees: {0}  VS  {1}".format(oldtree.GetEntries(),newtree.GetEntries()) )
    newtree.Write()

def update(runnumbers,channel_map_dict=None,TTree_dir="/Users/jeffburggraf/PycharmProjects/2NCorrrAnalysis/Data/TTrees/" ):
    """

    :param runnumbers:  runnumbers: the runnumbers of the root files. These file must be in TTRee_dir.
    :param channel_map_dict: See R1DC1990_conv()
    :param TTree_dir: Directory containing root files. The following naming convention is assumed: rXXXX.root
    :return: None
    """

    if isinstance(runnumbers, int):
        runnumbers = [runnumbers]

    for num in runnumbers:
        oldfilepath = TTree_dir + "r{}.root".format(num)
        newfilepath = TTree_dir + "newr{}.root".format(num)

        if not os.path.isfile(oldfilepath):
            assert False, "r{}.root does not exist".format(num)

        if not os.path.isfile(newfilepath):
            print "Creating TTree for run number {0} ... ".format(num)
            R1DC1990_conv(oldfilepath,channel_map_dict=channel_map_dict)
            print "Done.. TTree for run number {} up to date".format(num)
            continue

        oldfile = ROOT.TFile(oldfilepath)
        oldtree = oldfile.Get("R1DC")

        newfile = ROOT.TFile(newfilepath)
        newtree = newfile.Get("tree")

        oldentries = oldtree.GetEntries()
        newentries = newtree.GetEntries()

        assert oldentries != 0

        if newentries == oldentries:
            print "Run number {0} up to date with {1} pulses".format(num,
                                                                     newentries)

        if newentries < oldentries:
            print "Updating TTree for run # {0} from {1} entries to {2} entries ...".format(
                num, newentries, oldentries)
            R1DC1990_conv(oldfilepath, start_entry=newentries ,channel_map_dict=channel_map_dict)
            continue

def FF_conv(path_to_rootFile,start_entry=0,stop_entry=None,TDC_channel_map_dict=None,ADC_Cal=lambda x: 1*x+0, updating_file=False):
    """
    Restructures a daq2 r####.root TTree. Able to map selected R1DC1990[i][j]
    channels to a more user friendly leaf name. Resulting TTree MUCH faster and smaller.
    :param path_to_rootFile: Absolute pth to r####.root file.
    :param start_entry: Update
    :param stop_entry:
    :param channel_map_dict: Dictionary that maps R1DC1190[i][j]'s to new leaf names, e.g. R1DC1190[1][1]->D30T in the example below (See variable "dict" below. )
    :return:
    """

    def get_array_index(hit, rootVar):
        """TDC1190[100][129]/I<->TDC1190[hit][rootvar]/I   ... Used in R1DC1990()"""
        return hit * 129 + rootVar

    runnumber = int(
        re.search(r"r([0-9]+)\.root$", path_to_rootFile).groups(0)[0])

    # If root file not foud in given directory, check the download directory and move it to correct directory (symptom of laziness)
    if not os.path.isfile(path_to_rootFile):
        baseRootFileName=os.path.basename(path_to_rootFile)
        downloadsPath="/Users/jeffburggraf/Downloads/"
        if os.path.isfile(downloadsPath+baseRootFileName):

            copiedTreePath=os.path.dirname(path_to_rootFile)+"/"+baseRootFileName
            print downloadsPath+baseRootFileName,copiedTreePath

            copyfile(downloadsPath+baseRootFileName,copiedTreePath)

            assert os.path.isfile(copiedTreePath), "Coipying root file from download directory failed"

            FF_conv(copiedTreePath,start_entry,stop_entry,TDC_channel_map_dict)
            return None

        else:
            assert 0, "file {} cannot be found. ".format(path_to_rootFile)


    file=ROOT.TFile(path_to_rootFile)
    oldtree = file.Get("FF")
    if not isinstance(oldtree,ROOT.TTree):
        oldtree=file.Get("PAA")
        PAA=True
    else:
        PAA=False



    newpath = os.path.dirname(path_to_rootFile) + "/newr{0}.root".format(runnumber)

    if updating_file:
        assert os.path.isfile(newpath), "Cannot update {}, since it hasn't been created. "
        newfile = ROOT.TFile(newpath, "update")
        newtree = newfile.Get("tree")
        if oldtree.GetEntries() == newtree.GetEntries():
            return 1
        newtree.GetEntry(newtree.GetEntries())
        assert oldtree.GetEntries()>newtree.GetEntries()
        start_entry=newtree.GetEntries()
        print "Updating run number {0} from {1} entries to {2} entries...".format(runnumber,newtree.GetEntries(),oldtree.GetEntries())
    else:
        newfile=ROOT.TFile(newpath,"recreate")
        newtree=ROOT.TTree("tree","tree")


    if TDC_channel_map_dict==None:
        # Dictionary mapping R1DC1190[i][j] to new leaf name, e.g. R1DC1190[1][1]->D30T in the exaple below
        dict =  {(1, 1): "D30T", (1, 2): "D30B", (1, 3): "D90T", (1, 4): "D90B",
        (1, 5): "D150T", (1, 6): "D150B", (1, 7): "D210T",
        (1, 8): "D210B", (1, 9): "D120T", (1, 10): "D120B",
        (1, 11): "D270T", (1, 12): "D270B", (1, 13): "D330T",
        (1, 14): "D330B", (1,17):"HpGE_t",(1, 16): "trig"}
    else:
        dict=TDC_channel_map_dict

    TDC_addresses=[array.array("f",[0]) for key,value in sorted(dict.iteritems())]
    ADC_addresses=[array.array("f",[0]) for i in range(8)]

    TDC_names=[value for key,value in sorted(dict.iteritems())]
    ADC_names=["ADC{}".format(i) for i in range(8)]

    TDC_indicies=[get_array_index(key[0],key[1]) for key,value in sorted(dict.iteritems())]

    for TDC_name,TDC_address in zip(TDC_names,TDC_addresses):
        if updating_file:
            newtree.SetBranchAddress(TDC_name, TDC_address)
        else:
            newtree.Branch(TDC_name,TDC_address,TDC_name+"/F")

    for ADC_name, ADC_address in zip(ADC_names, ADC_addresses):
        if updating_file:
            newtree.SetBranchAddress(ADC_name, ADC_address)
        else:
            newtree.Branch(ADC_name, ADC_address, ADC_name + "/F")

    TDC_address_indicie=[(a,i) for a,i in zip(TDC_addresses,TDC_indicies)]

    ADC_address_indicie = [(a, i) for a, i in zip(ADC_addresses, range(8))]


    for i in xrange(start_entry,stop_entry if stop_entry else oldtree.GetEntries()):
        oldtree.GetEntry(i)

        if not PAA:
            for TDC_address,TDC_index in TDC_address_indicie:
                TDC_address[0]=oldtree.TDC1190[TDC_index]*0.1

        for ADC_address,ADC_index in ADC_address_indicie:
            ADC_address[0]=ADC_Cal(oldtree.PADC785N[ADC_index])

        newtree.Fill()


    if (oldtree.GetEntries()!=newtree.GetEntries() and stop_entry==None):
        warnings.warn( "Entry mismatch between trees: {0}  VS  {1}".format(oldtree.GetEntries(),newtree.GetEntries()) )
    newtree.Write()

class ENDF_crossection():
    """
    To access the properly formatted data for a given reaction, first plot the reaction, then click on "view evaluated data".
    """
    def __init__(self,path_to_evaluated_data_file, min=None, max = None):


        _file=open(path_to_evaluated_data_file)
        self.title = re.sub(";"," ",_file.readline())
        self.title=re.sub("&","#",self.title)
        self.title = re.sub(r"n<sub>0</sub>", r"n", self.title)
        self.title = re.sub(r"<sub>(.)</sub>",r"_{\1}",self.title)



        lines=[line.split()[0].split(",") for line in _file]

        self.xytuple = []

        for line in lines:
            try:
                result = tuple([float(line[0])/10.0**6, float(line[-1])])
            except:
                print "WARNING! Invalid line: {}\n".format(line)
                continue

            self.xytuple.append(result)

            if min is not None:
                if result[0] < min:
                    continue
            if max is not None:
                if result[0] > max:
                    break

        self.x,self.y = tuple(zip(*self.xytuple))

        self.TGraph=ROOT.TGraph(len(self.x),np.array(self.x,dtype=np.float64),np.array(self.y,dtype=np.float64))
        self.TGraph.SetTitle(self.title)
        self.TGraph.GetXaxis().SetTitle("[MeV]")
        self.TGraph.GetYaxis().SetTitle("#sigma [b]")

        self.min=self.x[0]
        self.max=self.x[-1]

    def Eval(self,x):
        assert self.min<=x<=self.max  , "Value {} is out of range! ".format(x)
        return self.TGraph.Eval(x)


class __negative_inf__():
    def __gt__(self, other):
        return False

    def __lt__(self, other):
        return True

    def __eq__(self, other):
        if isinstance(other, __negative_inf__):
            return True
        else:
            return False

class pad():
    __colors__=[("blue", 4), ("red", 2), ("black", 1),
     ("green", 3), ("cyan", 7), ("pink", 6),
     ("grey", 13), ("yellow", 5)]

    def __init__(self, pad):
        self.__Colors__ = [("blue", 4), ("red", 2), ("black", 1),
                       ("green", 3), ("cyan", 7), ("pink", 6),
                       ("grey", 13), ("yellow", 5)]

        self.graphColor_args=[]



        assert isinstance(pad, ROOT.TVirtualPad)
        self.title = None
        self.hist_titles = []
        self.hist_dict = {}
        self.hist_list=[]
        self.pad = pad
        self.max = __negative_inf__()
        self.min=10**8

        self.Legend=Legend()


        self.individualDrawOptions={} # Draw options  {Th1:option, ... }

        self.legend_list = []
        self.Legend_hist_dict={}
        self.LegendOptionDict={}
        self.main_hist = None

        self.drawables=[]
        self.draw_options=[]

        self.color_override_dict={}

        self.color_value_dict={}

        self.HasSetColors=0

    def add_drawable(self,object,options=""):

        self.drawables.append(object)
        self.draw_options.append(options)

    def cd(self):
        return self.pad.cd()

    def setTitle(self, title):
        self.title = title

    def get_hist(self,legend_entry):
        return self.Legend_hist_dict[legend_entry]

    def get_all_histos(self):
        return self.hist_dict.values()

    def add_graph_or_hist(self, hist_or_graph, legend_entry, LegendAddEntryOption="f", DrawOption="",LegTitle="",Override_Color=False):
        assert isinstance(hist_or_graph,(ROOT.TH1,ROOT.TGraph))
        _title=hist_or_graph.GetName()

        if DrawOption:
            self.individualDrawOptions[hist_or_graph]=DrawOption

        while _title  in self.hist_dict:
            _title=_title+"_"
            # hist_or_graph.SetTitle(_title)

        self.hist_titles.append(_title)


        self.hist_dict[_title] = hist_or_graph

        self.hist_list.append(hist_or_graph)

        if not self.main_hist:
            self.main_hist = self.hist_dict[self.hist_titles[
                0]]  # set main hist to first hist in list.

        if legend_entry:
            self.Legend.addentry(hist_or_graph, legend_entry, LegendAddEntryOption)
            if LegTitle:
                self.Legend.setTitle(LegTitle)

        if isinstance(hist_or_graph,ROOT.TH1):
            if hist_or_graph.GetLineColor()!=602:
                # print hist_or_graph.GetLineColor()
                Override_Color = hist_or_graph.GetLineColor()


        if Override_Color:
            self.color_override_dict[hist_or_graph]=Override_Color

    def __set_colors__(self):

        if self.HasSetColors:
            if len(self.hist_titles)==self.HasSetColors:
                return

        __ColorKeys__, __ColorValues__ = tuple(zip(*self.__Colors__))
        __ColorValues__ = list(__ColorValues__)

        line_styles = [1, 2, 3, 8]

        for index,hist_title in enumerate(self.hist_titles):
            hist=self.hist_dict[hist_title]

            if not isinstance(hist, (ROOT.TH1, ROOT.TGraph)): continue

            move_first_color_to_back=False


            if hist in self.color_override_dict:
                _color_=self.color_override_dict[hist]
            else:
                _color_=__ColorValues__[0]
                move_first_color_to_back=True


            if isinstance(hist,ROOT.TH1):
                if isinstance(hist,ROOT.TH1F):
                    if hist.GetLineColor()==602:
                        hist.SetLineColor(_color_)

                elif isinstance(hist,ROOT.TH2F):
                    hist.SetMarkerColor(_color_)
                    hist.SetLineColor(__ColorValues__[0])

            elif isinstance(hist, ROOT.TGraph):
                if hist.GetLineColor() == 1:
                    hist.SetLineColor(_color_)
            else:
                move_first_color_to_back=False

            self.color_value_dict[hist_title]=_color_  # use self.get_hist_color to retrieve this value


            if move_first_color_to_back==True:
                __ColorValues__ = __ColorValues__[1:] + [__ColorValues__[0]]

            if index >= len(__ColorValues__):
                if index >= 4 * len(__ColorValues__):
                    assert False, "Too many histograms."
                else:
                    if isinstance(hist,(ROOT.TH1F,ROOT.TGraph)):
                        hist.SetLineStyle(line_styles[index // len(__ColorValues__)])
                    else:
                        warnings.warn("Could not make unique color scheme for pad {}".format(self.title))

        self.HasSetColors=True

    # _pad.Draw(globalDrawOptions, sames, log, **legendOptions)
    def Draw(self, padWideDrawOptions="", sames=False, log=False,legendFontSize=0.05,**legendOptions):
        if not self.main_hist: return 0
        self.cd()


        self.__set_colors__()


        drawoptions = "".join(padWideDrawOptions.split())


        if log:
            ROOT.gPad.SetLogy()

        main_hist_draw_options=drawoptions
        if self.main_hist in self.individualDrawOptions:

            main_hist_draw_options+=self.individualDrawOptions[self.main_hist]
        if "a" not in main_hist_draw_options.lower() and isinstance(self.main_hist,ROOT.TGraph):
            main_hist_draw_options+="a"
        self.main_hist.Draw(main_hist_draw_options)

        # Set maximum and/or minimum (both for TGraph...)
        dont_set_range=False
        for hist_or_graph in self.hist_dict.values():
            if isinstance(hist_or_graph, ROOT.TH1F):
                histlist=HistToList(hist_or_graph)
                if log:
                    bin_values = filter(lambda x: x > 0, histlist.y)
                else:
                    bin_values=filter(lambda x:x!=0,histlist.y)

                if len(bin_values)==0:
                    continue

                _min = min(bin_values)
                _max= max(bin_values)

                if len(bin_values)>1:
                    if _max>self.max:
                        self.max=_max

                    if _min<self.min:
                        if (log and _min>0) or (log==False):
                            self.min=_min

                    if self.min==0:
                        warnings.warn("trouble determining minimum for {}. Setting min to 10^-6. ".format(hist_or_graph.GetTitle()))
                        self.min=10**-6

            elif isinstance(hist_or_graph,ROOT.TGraph):
                _yarray=np.frombuffer(hist_or_graph.GetY(), dtype=np.float64,count = hist_or_graph.GetN())
                _max = max(_yarray)
                _min=min(_yarray)
                if self.min==None:
                    self.min=_min

                if _max > self.max:
                    self.max = _max

                if _min < self.min:
                    self.min = _min
            else:
                dont_set_range=True


        if self.title:
            self.main_hist.SetTitle(self.title)

        if not (isinstance(self.max,__negative_inf__)):
            padding_range = abs(self.max - self.min)
        else:
            dont_set_range=True
            print "Not setting range for hist {}".format(hist_or_graph.GetTitle())


        if self.max==__negative_inf__:
            dont_set_range=True

        if dont_set_range==False:
            if log:
                if self.max>0 and self.min>0:

                    MMax=self.max**1.1/self.min**0.1
                    MMin = self.min ** 1.1 / self.max ** 0.1


                else:
                    dont_set_range=True
            else:
                MMax=(self.max)+(padding_range)*.1

                MMin=(self.min)-(padding_range)*.05

                if MMin<0 and self.min>0:
                    MMin=0

            if isinstance(hist_or_graph,(ROOT.TGraph,ROOT.TH1)) and self.max and dont_set_range==False:

                self.main_hist.SetMinimum(MMin)
                self.main_hist.SetMaximum(MMax)
            else:
                pass


        for hist_title in self.hist_titles[1:]:
            _hist = self.hist_dict[hist_title]
            _drawoptions = drawoptions

            if _hist in self.individualDrawOptions:
                _drawoptions+=self.individualDrawOptions[_hist]

            if isinstance(_hist, ROOT.TGraph):
                if "a" in _drawoptions:
                    _drawoptions = _drawoptions.replace("a", "")
            _hist.Draw(
                "same{}".format("s" if sames else "") + _drawoptions)

        # print "max:{0} = {1}".format(self.hist_titles, self.max, )
        self.Legend.Draw(fontsize = legendFontSize,**legendOptions)
            # self.LG = Legend(self.legend_list, **legendOptions)

        for drawable,option in zip(self.drawables,self.draw_options):
            drawable.Draw(option)


        self.pad.Update()
        self.pad.Modified()

__thesis_plot_refs__ = []
def thesis_plot(hisos_and_or_TGraphsList, LegendEntries=[],TGraphOptions="*",TH1options="hist",draw=True,legfontsize=0.03,ndigits = 2):
    ROOT.TGaxis.SetMaxDigits(ndigits)

    if draw:
        c1=ROOT.TCanvas("canvas","Canvas Title",1000, 600)
        c1.SetFixedAspectRatio(0)
        __thesis_plot_refs__.append(c1)

    if LegendEntries:
        assert len(hisos_and_or_TGraphsList)==len(LegendEntries)

    try:
        ROOT.gPad.SetLeftMargin(0.1)
        ROOT.gPad.SetBottomMargin(0.15)
    except:
        __thesis_plot_refs__.append(ROOT.TCanvas())
        ROOT.gPad.SetLeftMargin(0.1)
        ROOT.gPad.SetBottomMargin(0.15)

    leg =ROOT.TLegend(0.8,0.8,0.95,0.95)
    leg.SetTextSize(legfontsize)
    leg.SetLineColor(ROOT.kWhite)
    __thesis_plot_refs__.append(leg)

    for i,drawable in enumerate(hisos_and_or_TGraphsList):
        xaxis = drawable.GetXaxis()
        yaxis = drawable.GetYaxis()
        xaxis.CenterTitle();
        yaxis.CenterTitle();
        xaxis.SetTitleOffset(1.5)
        yaxis.SetTitleOffset(1.3)
        xaxis.SetTitleSize(0.04)
        yaxis.SetTitleSize(0.04)
        xaxis.SetLabelSize(0.03)
        yaxis.SetLabelSize(0.03)

        if LegendEntries:
            leg.AddEntry(drawable,LegendEntries[i],"LP")
        if draw:
            if isinstance(drawable,ROOT.TGraph):
                if i:
                    drawable.Draw("same "+TGraphOptions)
                else:
                    drawable.Draw("A" + TGraphOptions)
            elif isinstance(drawable,ROOT.TH1):
                if i:
                    drawable.Draw("same "+TH1options)
                else:
                    drawable.Draw(TH1options)
            else:
                assert 0

    if LegendEntries:
        leg.Draw()



class multipad():
    num_calls=1
    """
    Easily create and divide a canvas and place multiple histograms within.
    Example:
        padnames=[30,90,120,150,225,270] # Name of all the pad's.
        m=multipad(padnames)

        padtitles=["{}_title".format(name) for name in padnames ] # Titles that will appear above histos for each pad
        m.set_pad_titles(padtitles)  # set pad titles with list. Use same order as padnames.

        for padname in padnames:
            for j in range(4):
                _name="{0}{1}".format(padname,j) # histogram name
                _hist=ROOT.TH1F(_name,_name,20,-10,10)  # Reference to hist
                _hist.FillRandom("gaus",(j+1)*500)

                # add the histogram to the appropriate pad ,and give our histogram an entry in the legend
                m.add_hist(_hist, padname,legend_entry=str(j)+"_Legend")

        m.Draw(drawoptions="E hist",sames=True,fontsize=0.02)

    """

    def __init__(self, Pad_titles_list, LegendTitleList=None ,CanvasTitle="auto"):

        if not isinstance(Pad_titles_list,list):
            Pad_titles_list=[Pad_titles_list]

        self.LegendTitleList = LegendTitleList
        if self.LegendTitleList:
            if not isinstance(self.LegendTitleList, list):
                self.LegendTitleList = [self.LegendTitleList]
            assert len(self.LegendTitleList)==len(Pad_titles_list) , 'Number of legend titles must equal number of pad titles. '

        assert isinstance(Pad_titles_list, list)
        self.global_hist_max=0

        Pad_titles_list = map(str, Pad_titles_list)

        n = len(Pad_titles_list)
        if n < 4:
            divide = (n,)
        elif n == 4:
            divide = (2, 2)
        elif n == 5 or n == 6:
            divide = (3, 2)
        elif n == 7 or n == 8:
            divide = (4, 2)
        elif n == 9:
            divide = (3, 3)
        elif 10 <= n <= 12:
            divide = (4, 3)
        elif 13 <= n <= 16:
            divide = (4, 4)
        else:
            assert False, "Too many pads to fit on a single canvas! Max is 16 pads. "

        y = 450 if len(divide) == 1 else 900

        if CanvasTitle=="auto":
            CanvasTitle="c"+str(multipad.num_calls)
        self.canvas = ROOT.TCanvas(CanvasTitle, CanvasTitle, 1220,
                                   y)
        self.canvas.Divide(*divide)

        self.padTitles = Pad_titles_list
        self.pad_dict = OrderedDict()

        for index, title in enumerate(self.padTitles):
            _pad = self.canvas.cd(index + 1)
            _pad.SetTitle(title)
            self.pad_dict[title] = pad(_pad)

        multipad.num_calls+=1

    def update(self):
        for _pad in self.get_pads():
            _pad.pad.Update()
            _pad.pad.Modified()

    def get_color(self,hist_title_or_TH1):

        failed_to_find=0

        if isinstance(hist_title_or_TH1,ROOT.TH1):
            return hist_title_or_TH1.GetLineColor()
        elif isinstance(hist_title_or_TH1,(ROOT.TH1,ROOT.TGraph)):
            hist_title_or_TH1=hist_title_or_TH1.GetTitle()

        hist_title=hist_title_or_TH1

        for _pad in self.get_pads():
            assert isinstance(_pad,pad)
            _pad.__set_colors__()

            if failed_to_find:
                break

            for _hist in _pad.get_all_histos():
                if _hist.GetTitle()==str(hist_title):
                    if len(_pad.color_value_dict.values()) == 0:
                        failed_to_find = 2
                        break


                    if hist_title in _pad.color_value_dict:
                        return _pad.color_value_dict[hist_title]
                    else:
                        failed_to_find=1
            if failed_to_find:
                break
        else:
            failed_to_find=1

        if failed_to_find:
            if failed_to_find==2:
                warnings.warn(
                    "Could not find line color for {}. Be sure to do plot.Draw() first!".format(hist_title))
            else:
                warnings.warn("Could not find line color for {} .".format(hist_title))
            return 0


    def get_pads(self):
        return [self.pad_dict[t] for t in self.padTitles]

    def set_pad_margins(self,left=0.1,right=0.1):
        for pad in self.get_pads():
             pad.pad.SetLeftMargin(left);
             pad.pad.SetRightMargin(right);
             pad.pad.Modified()
             pad.pad.Update()
        self.canvas.Update()

    def SavePicture(self,name,dir=""):
        name=str(name)
        if dir[-1]!="/":
            dir+="/"
        if name[-3:]!=".ps":
            name+=".ps"

        self.canvas.Update()
        self.canvas.SaveAs(dir+name)

    def saveToROOTFile(self,file_name,dir):
        if file_name[-5:]!=".root":
            file_name+=".root"
        if dir[-1]!="/":dir+="/"
        newfile=ROOT.TFile(dir+file_name,"recreate")
        newfile.cd()

        for _hist in self.get_all_histos():
            _hist.Write()
        newfile.Close()

    def set_colors(self):
        for _pad in self.get_pads():
            assert isinstance(_pad,pad)
            _pad.__set_colors__()

    def set_pad_titles(self, title_list):
        assert len(title_list) == len(self.padTitles)
        for title, pad in zip(title_list, self.get_pads()):
            pad.setTitle(title)

    def add_drawables(self,drawable, to_pad_with_title,draw_options=""):
        assert hasattr(drawable,"Draw"), "Object {} does not have Draw() method. ".format(drawable)
        to_pad_with_title = str(to_pad_with_title)
        assert to_pad_with_title in self.padTitles, " pad titled {} does not exist".format(
            to_pad_with_title)
        _pad = self.pad_dict[to_pad_with_title]
        assert isinstance(_pad, pad)

        _pad.add_drawable(drawable,options=draw_options)


    def add_hist(self, hist_or_graph, to_pad_with_title, legend_entry=None,legend_entryOption="f",DrawOption="",override_title=False,Override_Color=False):
        """
        :param hist:  ROOT TH1 or similar...
        :param to_pad_with_title: pad title as set in multipad constructor
        :param legend_entry: ...
        :param legend_entryOption: Option for TLegend.AddEntry.
        :param DrawOption: Draw option used for this histogram only.
        :return: ROOT.TPad that hist was added to.
        """
        if isinstance(hist_or_graph,ROOT.TH1):
            _max=hist_or_graph.GetBinContent(hist_or_graph.GetMaximumBin())
            if _max>self.global_hist_max:
                self.global_hist_max=_max

        to_pad_with_title = str(to_pad_with_title)
        assert to_pad_with_title in self.padTitles, " pad titled {} does not exist".format(
            to_pad_with_title)
        _pad = self.pad_dict[to_pad_with_title]
        assert isinstance(_pad, pad)

        if self.LegendTitleList:
            LegTitle = self.LegendTitleList[
                self.padTitles.index(to_pad_with_title)]
        else:
            LegTitle=""


        _pad.add_graph_or_hist(hist_or_graph, legend_entry, legend_entryOption, DrawOption,LegTitle,Override_Color=Override_Color)

        if override_title:
            _pad.main_hist.SetTitle(str(override_title))
        return _pad.pad

    def Draw(self, globalDrawOptions="", sames=False,log=False, use_global_Y_max=0,gridlines=False, legendFontSize=0.05,**legendOptions):

        if hasattr(self, "__thesisCalled__"):
            legendFontSize=0.03

        for pad_title in self.padTitles:
            _pad = self.pad_dict[pad_title]
            assert isinstance(_pad, pad)
            _pad.Draw(globalDrawOptions, sames, log,legendFontSize, **legendOptions)

        if self.global_hist_max and use_global_Y_max:
            for _hist in self.get_all_histos():
                _hist.SetMaximum(self.global_hist_max * 1.1)

        for _pad in self.get_pads():
            _pad.pad.Modified()
            if gridlines:
                _pad.pad.SetGridx()
                _pad.pad.SetGridy()

        # if log:  # change min/max for log plot
        #     for _pad in self.get_pads():
        #         assert isinstance(_pad,pad)
        #         print _pad.max,_pad.min
        #
        #         _pad.main_hist.SetMaximum(_pad.max**1.1/_pad.min**0.1)



        self.canvas.Update()

        if hasattr(self,"__thesisCalled__"):
            self.__thesis_post_Draw__()

    def set_global_min_max(self,Min=None,Max=None):
            for _pad in self.get_pads():
                assert isinstance(_pad,pad)
                if isinstance(_pad.main_hist,(ROOT.TH1,ROOT.TH2)):
                    if Min:
                        _pad.main_hist.SetMinimum(Min)
                    if Max:
                        _pad.main_hist.SetMaximum(Max)

    def normalize_all_histos(self):
        for _h in self.get_all_histos():
            if isinstance(_h,(ROOT.TH1,ROOT.TH2)):
                normalize_hist(_h)

    def thesis(self):
        self.__thesisCalled__=True
        ROOT.c1.SetFixedAspectRatio(0)
        ROOT.c1.SetWindowSize(1000, 600)
        for _pad in self.get_pads():
            histo = _pad.main_hist
            xaxis = histo.GetXaxis()
            yaxis = histo.GetYaxis()
            xaxis.CenterTitle();
            yaxis.CenterTitle();
            xaxis.SetTitleOffset(1.5)
            yaxis.SetTitleOffset(1.5)
            xaxis.SetTitleSize(0.04)
            yaxis.SetTitleSize(0.04)
            xaxis.SetLabelSize(0.03)
            yaxis.SetLabelSize(0.03)
            _pad.pad.SetLeftMargin(0.12)
            _pad.pad.SetBottomMargin(0.2)

            _pad.pad.Update()


    def __thesis_post_Draw__(self):
        for _pad in self.get_pads():
            _pad.Legend.leg

    def get_all_histos(self):
        pads=self.get_pads()
        histos=[]
        for pad in pads:
            histos+=pad.get_all_histos()
        return histos

    def get_histo(self, pad_title, hist_title=None):
        pad_title=str(pad_title)
        assert pad_title in self.padTitles, "Pad titled {0} does not exist. \n Titles: \n {1}".format(
            pad_title,self.padTitles)
        _pad = self.pad_dict[pad_title]
        assert isinstance(_pad, pad)

        if hist_title:
            assert hist_title in _pad.hist_titles, "hist titled {0} in pad {1} does not exist. Titles: \n {2}".format(
                hist_title,pad_title,_pad.hist_dict.keys())
            return _pad.hist_dict[hist_title]

        else:
            return [_pad.hist_dict[title] for title in
                    _pad.hist_titles]

    def get_main_histos(self):
        return [_pad.main_hist for _pad in self.get_pads()]

    def getLegendsList(self):
        return [_pad.Legend for _pad in self.pad_dict.values()]

def cut_allNoneZero(*strings_or_list):
    vars=[]
    for thing in strings_or_list:
        if hasattr(thing, "__iter__"):
            for var in thing:
                vars.append(str(var))
        else:
            vars.append(str(thing))

    return "!=0 && ".join(vars)+"!=0"

def cut_rangeAND(Range_tuple,*ANDexpressions):
    cuts=[]
    expressions=flatten_and_concatenate(ANDexpressions)
    for thing in expressions:
        thing=str(thing)
        cut1="("+thing+")"+">"+str(Range_tuple[0])
        cut2=" && ("+thing+")"+"<"+str(Range_tuple[1])
        cuts.append(cut1+cut2)
    return " && ".join(cuts)

def cut_rangeOR(Range_tuple,*ORexpressions):
    cuts=[]
    expressions=flatten_and_concatenate(ORexpressions)
    for thing in expressions:
        thing=str(thing)
        cut1="("+thing+")"+">"+str(Range_tuple[0])
        cut2=" && ("+thing+")"+"<"+str(Range_tuple[1])
        cuts.append("("+cut1+cut2+")")
    return " || ".join(cuts)

def cut_range_outer(range_tuple,*expressions):
    range_tuple=sorted(range_tuple)[::-1]
    cuts = []
    expressions = flatten_and_concatenate(expressions)
    for thing in expressions:
        cut1="({0})>{1}".format(thing,range_tuple[0])
        cut2="({0})<{1}".format(thing,range_tuple[-1])
        cuts.append("({0} || {1})".format(cut1,cut2))
    return "&&".join(cuts)

def cut_AND(*List):
    List=flatten_and_concatenate(List)
    return "("+"&& ".join(List)+")"

def cut_OR(*List):
    List = flatten_and_concatenate(List)
    return "("+"||".join(List)+")"



def getToF(degree):
    return "(D{0}T+D{0}B)/2-trig".format(degree)

def cut_TrigTopBot(degree):
    return "(trig && D{0}T && D{0}B)".format(degree)

def normalize_graph(tgraphs):

    for graph in tgraphs:
        assert isinstance(graph,ROOT.TGraph)
        x = np.array(graph.GetX(), dtype=np.float32)
        y = np.array(graph.GetY(), dtype=np.float32)
        errx = np.array(graph.GetEX(), dtype=np.float32)
        erry = np.array(graph.GetEY(), dtype=np.float32)
        norm = sum(y)

        y/=norm
        erry/=norm

        for i in range(len(x)):
            graph.SetPoint(i,x[i],y[i])
            if isinstance(graph,ROOT.TGraphErrors):
                graph.SetPointError(i,errx[i],erry[i])

def normalize_hist(hist_OR_listOfHists,norm_or_TH1F=1, UseBinWidthForDx=True, SquareIntegral=False):

    if not hasattr(hist_OR_listOfHists,"__iter__"):
        hist_OR_listOfHists = [hist_OR_listOfHists]

    for hist in hist_OR_listOfHists:
        if SquareIntegral:
            hlist = HistToList(hist)
            y = np.array(hlist.y)
            mag = sum(y*y)*hist.GetBinWidth(1)
            hist.Sumw2()
            hist.Scale(1.0 / np.sqrt(mag))
        else:
            if isinstance(norm_or_TH1F,ROOT.TH1):
                norm_or_TH1F=norm_or_TH1F.Integral()

            mag=float(hist.Integral("width" if UseBinWidthForDx else ""))

            if mag==0:
                warnings.warn("Histogram {} not normailized because integral was zero.".format(hist.GetTitle()))
                continue

            hist.Scale(float(norm_or_TH1F)/mag)
            hist.Sumw2()

    return hist_OR_listOfHists


def GetMostRecentHist():
    return ROOT.gPad.GetListOfPrimitives()[-1]

class bin_info():
    def __init__(self,hist):
        self.hist=hist
        self.axis=[self.hist.GetXaxis()]

        if isinstance(hist,ROOT.TH2):
            self.axis.append(self.hist.GetYaxis())

        if isinstance(hist, ROOT.TH3):
            self.axis.append(self.hist.GetZaxis())

        self.widths=map(lambda x:x.GetBinWidth(1),self.axis)
        self.nbins = map(lambda x: x.GetNbins(), self.axis)
        self.mins= map(lambda x:x.GetXmin(),self.axis)

    def get_bin_number_tuple(self,*cords):
        # assert len(cords)==len(self.widths) , "must supply {} coordinates".format(len(self.widths))
        N=len(cords)

        result=[int((xi-x0)/w+1) for x0,xi,w in zip(self.mins[:N],cords,self.widths[:N]) ]

        for index,(binnum,maxbinNum) in enumerate(zip(result,self.nbins[:N])):
            if binnum>maxbinNum:
                result[index]=maxbinNum+1
            if binnum<0:
                result[index]=0

        return tuple(result)

    def get_bin_number_scalar(self, *cords):
        return self.hist.GetBin(*self.get_bin_number_tuple(*cords))


# class TH1F(ROOT.TH1F):
#     def __init__(self, name="", title="", nbins=1, xmin=0, xmax=1,TH1F=None, *args):
#         ROOT.TH1F.__init__(self, name, title, nbins, xmin, xmax, *args)
#
#
#         self.binwidth=self.GetXaxis().GetBinWidth(1)
#         self.xmin=xmin;self.xmax=xmax;self.nbins=nbins
#
#         self.underflow=[];self.overflow=[]
#
#         self.x=np.arange(xmin+self.binwidth/2.0,xmax+self.binwidth/2.0,float(xmax-xmin)/nbins)
#
#         self.bin_content=np.array([0.0]*(nbins+2),dtype=np.float64)
#
#     def getTGraph(self):
#         self.graph=ROOT.TGraph(self.nbins, self.x,  self.bin_content[1:-1])
#         return self.graph
#
#     def Eval(self, x):
#         return self.Interpolate(x)
#
#     def Fill(self,x,w=1,*args):
#         binnum=int(float(x)/self.binwidth)+1
#         if binnum>self.nbins:binnum=self.nbins+1
#         elif binnum<0:binnum=0
#
#         self.bin_content[binnum]+=w
#
#         ROOT.TH1F.Fill(self,x,w,*args)
#
#     def Clone(self, *args):
#         return ROOT.TH1F.Clone(self, *args)
#
#     def Copy(self, *args):
#         return ROOT.TH1F.Copy(self, *args)
#
#     def Divide(self, *args):
#         return ROOT.TH1F.Divide(self, *args)
#
#     def Draw(self, *args):
#         return ROOT.TH1F.Draw(self, *args)
#
#     def FillRandom(self, *args):
#         return ROOT.TH1F.FillRandom(self, *args)
#
#     def FindBin(self, *args):
#         return ROOT.TH1F.FindBin(self, *args)
#
#     def Fit(self, *args):
#         return ROOT.TH1F.Fit(self, *args)
#
#     def GetBin(self, *args):
#         return ROOT.TH1F.GetBin(self, *args)
#
#     def GetBinCenter(self, binnum):
#         return ROOT.TH1F.GetBinCenter(self, binnum)
#
#     def GetBinContent(self, binnum):
#         return ROOT.TH1F.GetBinContent(self,binnum)
#
#     def GetBinError(self, binnum):
#         return ROOT.TH1F.GetBinError(self,binnum)
#
#     def GetBinLowEdge(self, binnum):
#         return ROOT.TH1F.GetBinLowEdge(self, binnum)
#
#     def GetBinWidth(self, binnum=1):
#         return ROOT.TH1F.GetBinWidth(self,binnum)
#
#     def GetBinWithContent(self, c, binx, firstx=0, lastx=0, maxdiff=0):
#         """Compute first binx in the range [firstx,lastx] for which diff = abs(bin_content-c) <= maxdiff."""
#         return ROOT.TH1F.GetBinWithContent(self, c, binx, firstx, lastx, maxdiff)
#
#     def GetEntries(self):
#         return ROOT.TH1F.GetEntries(self)
#
#     def GetMaximumBin(self,*args):
#         return ROOT.TH1F.GetMaximumBin(self,*args)
#     def GetMean(self,*args):
#         return ROOT.TH1F.GetMean(self,*args)
#     def GetMeanError(self,*args):
#         return ROOT.TH1F.GetMeanError(self,*args)
#     def GetMinimumBin(self,*args):
#         return ROOT.TH1F.GetMinimumBin(self,*args)
#     def GetName(self,*args):
#         return ROOT.TH1F.GetName(self,*args)
#     def GetNbinsX(self,*args):
#         return ROOT.TH1F.GetNbinsX(self,*args)
#     def GetRMS(self,*args):
#         return ROOT.TH1F.GetRMS(self,*args)
#     def GetStats(self,*args):
#         return ROOT.TH1F.GetStats(self,*args)
#     def GetStdDev(self,*args):
#         return ROOT.TH1F.GetStdDev(self,*args)
#     def GetStdDevError(self,*args):
#         return ROOT.TH1F.GetStdDevError(self,*args)
#     def GetTitle(self,*args):
#         return ROOT.TH1F.GetTitle(self,*args)
#     def Integral(self,*args):
#         return ROOT.TH1F.Integral(self,*args)
#     def Multiply(self,*args):
#         return ROOT.TH1F.Multiply(self,*args)
#     def Rebin(self,*args):
#         return ROOT.TH1F.Rebin(self,*args)
#     def RebinAxis(self,*args):
#         return ROOT.TH1F.RebinAxis(self,*args)
#     def RebinX(self,*args):
#         return ROOT.TH1F.RebinX(self,*args)
#     def Scale(self,*args):
#         return ROOT.TH1F.Scale(self,*args)
#     def SetAxisColor(self,*args):
#         return ROOT.TH1F.SetAxisColor(self,*args)
#     def SetAxisRange(self,*args):
#         return ROOT.TH1F.SetAxisRange(self,*args)
#     def SetBarOffset(self,*args):
#         return ROOT.TH1F.SetBarOffset(self,*args)
#     def SetBarWidth(self,*args):
#         return ROOT.TH1F.SetBarWidth(self,*args)
#     def SetBinContent(self,*args):
#         return ROOT.TH1F.SetBinContent(self,*args)
#     def SetBinError(self,*args):
#         return ROOT.TH1F.SetBinError(self,*args)
#     def SetBinErrorOption(self,*args):
#         return ROOT.TH1F.SetBinErrorOption(self,*args)
#     def SetDrawOption(self,*args):
#         return ROOT.TH1F.SetDrawOption(self,*args)
#     def SetFillAttributes(self,*args):
#         return ROOT.TH1F.SetFillAttributes(self,*args)
#     def SetFillColor(self,*args):
#         return ROOT.TH1F.SetFillColor(self,*args)
#     def SetFillColorAlpha(self,*args):
#         return ROOT.TH1F.SetFillColorAlpha(self,*args)
#     def SetFillStyle(self,*args):
#         return ROOT.TH1F.SetFillStyle(self,*args)
#     def SetLabelColor(self,*args):
#         return ROOT.TH1F.SetLabelColor(self,*args)
#     def SetLabelFont(self,*args):
#         return ROOT.TH1F.SetLabelFont(self,*args)
#     def SetLabelOffset(self,*args):
#         return ROOT.TH1F.SetLabelOffset(self,*args)
#     def SetLabelSize(self,*args):
#         return ROOT.TH1F.SetLabelSize(self,*args)
#     def SetLineColor(self,*args):
#         return ROOT.TH1F.SetLineColor(self,*args)
#     def SetLineColorAlpha(self,*args):
#         return ROOT.TH1F.SetLineColorAlpha(self,*args)
#     def SetLineStyle(self,*args):
#         return ROOT.TH1F.SetLineStyle(self,*args)
#     def SetLineWidth(self,*args):
#         return ROOT.TH1F.SetLineWidth(self,*args)
#     def SetMarkerAttributes(self,*args):
#         return ROOT.TH1F.SetMarkerAttributes(self,*args)
#     def SetMarkerColor(self,*args):
#         return ROOT.TH1F.SetMarkerColor(self,*args)
#     def SetMarkerColorAlpha(self,*args):
#         return ROOT.TH1F.SetMarkerColorAlpha(self,*args)
#     def SetMarkerSize(self,*args):
#         return ROOT.TH1F.SetMarkerSize(self,*args)
#     def SetMarkerStyle(self,*args):
#         return ROOT.TH1F.SetMarkerStyle(self,*args)
#     def SetMaximum(self,*args):
#         return ROOT.TH1F.SetMaximum(self,*args)
#     def SetMinimum(self,*args):
#         return ROOT.TH1F.SetMinimum(self,*args)
#     def SetName(self,*args):
#         return ROOT.TH1F.SetName(self,*args)
#     def SetNdivisions(self,*args):
#         return ROOT.TH1F.SetNdivisions(self,*args)
#     def SetStats(self,*args):
#         return ROOT.TH1F.SetStats(self,*args)
#     def SetTickLength(self,*args):
#         return ROOT.TH1F.SetTickLength(self,*args)
#     def SetTitle(self,*args):
#         return ROOT.TH1F.SetTitle(self,*args)
#     def SetTitleFont(self,*args):
#         return ROOT.TH1F.SetTitleFont(self,*args)
#     def SetTitleOffset(self,*args):
#         return ROOT.TH1F.SetTitleOffset(self,*args)
#     def SetTitleSize(self,*args):
#         return ROOT.TH1F.SetTitleSize(self,*args)
#     def SetXTitle(self,*args):
#         return ROOT.TH1F.SetXTitle(self,*args)
#     def SetYTitle(self,*args):
#         return ROOT.TH1F.SetYTitle(self,*args)
#     def ShowPeaks(self,*args):
#         return ROOT.TH1F.ShowPeaks(self,*args)
#     def Smooth(self,*args):
#         return ROOT.TH1F.Smooth(self,*args)


def event_rate_histogram(tree,delta_t=1,pulse_rate_Hz=240,cut=""):
    events_per_fill=pulse_rate_Hz*delta_t
    number_of_fills=tree.GetEntries()//events_per_fill

    fills=np.zeros(number_of_fills,dtype=np.intc)
    i=0

    while i<len(fills):
        num_hits=tree.Project("","",cut,"",events_per_fill,i*events_per_fill)
        fills[i]=num_hits
        i+=1
    _min=min(fills)-np.std(fills)
    _max=max(fills)+np.std(fills)
    nbins=int(float(_max-_min)/(0.25*np.std(fills)))
    h=ROOT.TH1F("h","h",nbins,_min,_max)
    [h.Fill(fills[i]) for i in xrange(len(fills))]
    return h

def time_dependant_hist(tree, delta_t_per_bin_minutes, cut="",
                        hist_title=None, pulse_rate_Hz=240, max_entry=None,n_pulses_overwrite=None,stats=0):
    assert isinstance(tree, ROOT.TTree)

    if max_entry == None:
        max_entry = tree.GetEntries() if n_pulses_overwrite==None else n_pulses_overwrite


    nbins = int(
        (max_entry) / float(pulse_rate_Hz) / (delta_t_per_bin_minutes * 60))
    if nbins==0:nbins=20
    total_time = max_entry/pulse_rate_Hz/60.
    bin_width_pulses = int(tree.GetEntries() / float(nbins))

    if hist_title == None:
        if not hasattr(time_dependant_hist, "hist_title_index"):
            time_dependant_hist.hist_title_index = 0
        else:
            time_dependant_hist.hist_title_index += 1

        hist_title = "___" + str(time_dependant_hist.hist_title_index)

    hist_title = str(hist_title)
    _hist = ROOT.TH1F(hist_title, hist_title, nbins, 0, total_time)

    bin_centers = HistToList(_hist).x

    pulse_ranges = [(i * bin_width_pulses, (i + 1) * bin_width_pulses) for i
                    in range(nbins)]

    for bin_center, (min_evt, max_evt) in zip(bin_centers, pulse_ranges):
        num_evts = tree.Project("", "", cut, "", max_evt - min_evt, min_evt)
        _hist.Fill(bin_center, num_evts)

    for i in range(1, nbins + 1):
        _hist.SetBinError(i, 0)
    _hist.SetStats(stats)

    return _hist









# bottom
"""

if __name__=="__main__":
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


"""
