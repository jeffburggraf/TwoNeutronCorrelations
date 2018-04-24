import ROOT
import numpy as np
import numbers
from collections import OrderedDict
from itertools import product
import time
from itertools import izip
import inspect
import itertools
from mytools import round as mtround
import warnings

def apply_operator_along_axis(arr, b, princ_axis=None, operator = lambda a,b:a*b):
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr)
    if not isinstance(b, np.ndarray):
        b = np.array(b)

    assert len(inspect.getargspec(operator).args)==2, "Operator must be a binary function which takes exactly two arguments."

    if princ_axis is None:
        assert arr.shape == b.shape,\
            "\nIf no princ_axis argument is supplied, then arrays must have the same shape. "
        return arr*b

    _princ_axis_values = ["x","y","z","xy","xz","yz"]

    assert isinstance(princ_axis, str), \
        "\n'princ_axis' must be a string. The valid values for princ_axis are: {}".format("'"+"', '".join(_princ_axis_values)+"'")

    princ_axis = princ_axis.lower()

    assert  princ_axis in ["x","y","z","xy","xz","yz"],\
        "\n Invalid 'princ_axis'. The valid values for princ_axis are: {}".format("'"+"', '".join(_princ_axis_values)+"'")

    princ_str = princ_axis
    princ_axis_indicies = ["xyz".index(i) for i in princ_axis]

    assert len(princ_axis_indicies)==b.ndim, \
        "\nToo many axis ({1} provided) to use with a {0}-dim array".format(b.ndim,len(princ_axis_indicies))

    assert max(princ_axis_indicies)<arr.ndim, "The array being multiplied by does not have a {0} axis!".format("xyz"[max(princ_axis_indicies)])

    assert all([arr.shape[i_a]==b.shape[i_other] for i_other,i_a in enumerate(princ_axis_indicies)]), \
        "\nThe {0} axis of array a ({0} axis shape is {1}) does not match the shape of b (b's shape is {2})!".format("-".join(princ_str),tuple(np.array(arr.shape)[princ_axis_indicies]), b.shape)

    broadcaster = [(slice(None) if i in princ_axis_indicies else np.newaxis) for i in range(arr.ndim)]

    broadcasted_b = np.broadcast_to(b[tuple(broadcaster)],arr.shape)
    result = map(lambda x: operator(*x), zip(arr.ravel(),broadcasted_b.ravel()))
    result = np.array(result).reshape(arr.shape)

    return result


def __parse_multidim_args_befor_hist__(ndim, mins, maxs, nbinss, binwidths, binarrays):
    """ Use to prepare the args and check for errors before passing them to the <hist> class. """
    for arg, arg_name in zip([mins, maxs, nbinss, binwidths],["mins", "maxs", "nbinss","binwidths"]):
        if arg is not None:
            assert  ((hasattr(arg, "__iter__") and len(arg) == ndim)) or  isinstance(arg,numbers.Number),\
            "\nArgument <{0}> must be either a single number (assumed to be the same for all axis), or an aray of len {1}\nMaybe you meant to use TH{2}F?".format(arg_name,ndim,len(arg))

            if not hasattr(arg, "__iter__"):
                locals()[arg_name] =  np.array([0]*ndim, dtype=np.float32)
                exec("{0} = np.array([{0}]*ndim, dtype=np.float32)".format(arg_name))
            else:
                exec("{0} = np.array({0}, dtype=np.float32)".format(arg_name))

    if binarrays is not None:
        assert hasattr(binarrays,"__iter__"), "binarrays must be iterable!"
        if all([hasattr(i,"__iter__") for i in binarrays]):
            assert len(binarrays)==ndim, "Invalid number of bin arrays in <bionarrays> argument! <binarrays> may either be a list of {0} arrays," \
                                         " or a single array which will be used for each of the {0} axis.".format(ndim)
        elif  all([isinstance(i,numbers.Number) for i in binarrays]):
            binarrays = [binarrays[:]]*ndim
        else:
            assert False, "Incorrect xbins specification! "

    if ((mins is not None) and (maxs is not None) and (nbinss is None) and (binwidths is not None)):
        nbinss =  np.array((maxs-mins)/binwidths,dtype=np.int)
        for i, _nbin_ in enumerate(nbinss):
            assert _nbin_>0, "The bin width, {0}, is too big for its corresponding range, [{1:1.2f}, {2:1.2f}]".format(binwidths[i], mins[i],maxs[i])

    return (mins, maxs, nbinss, binarrays)


def parrseHistArgs( mins, maxs, nbins, binarrays, name, title):
    """used by the hist class to parse args"""
    args = OrderedDict()

    if name is None:
        name = "Histo{0}".format(hist.name_index)
        hist.name_index += 1
    elif name in hist.names_used:
        hist.names_used[name] += 1
        name += str(hist.names_used[name]-1)
    else:
        hist.names_used[name] = 1

    if title == None:
        title = ""

    args["name"] = name
    args["title"] = title

    if mins is not None and maxs is not None:

        assert binarrays is None, "Cannot specify min/max along with an array of bin left edges.  "
        assert nbins is not None, "must specify number of bins!"

        assert hasattr(mins, "__iter__") and hasattr(maxs, "__iter__") and hasattr(nbins,
                                                                                   "__iter__"), "Args must be iterable"

        assert len(mins) == len(maxs) == len(nbins), "All args must be of same dimension!"
        ndim = len(mins)

        for i, (_nbins_, _min_, _max_) in enumerate(zip(nbins, mins, maxs)):
            assert all([isinstance(inst, numbers.Number) for inst in
                        [_min_, _max_, _nbins_]]), "Bin spec args must contain only numbers.  "
            assert _nbins_ > 0, "nbins is not greater than zero!"
            args["nbins{}".format(i)] = int(_nbins_)
            args["min{}".format(i)] = _min_
            args["max{}".format(i)] = _max_

    elif mins is None and maxs is None and nbins is None and binarrays is not None:
        ndim = len(binarrays)
        assert ndim in [1, 2, 3], "Wrong number of dimensions in xbins arrays"

        for i, xbins in enumerate(binarrays):
            args["nbins{}".format(i)] = len(xbins) - 1
            args["bins{}".format(i)] = np.array(xbins, dtype=np.float32)
            __bins__ = args["bins{}".format(i)]
            assert len(__bins__.shape)==1, "Arrays used to specify bins must be one-dimensional. "

    else:
        assert False, "Arguments provided are suported"

    return args,ndim

class hist(object):
    tcanvii_refs = []
    name_index = 0
    names_used = {} # dictionary of names and number of times used, e.g. {"some_name":2,...}

    def __init__(self, mins=None, maxs=None, nbins=None, binarrays=None, name=None, title=None):

        self.__self_args__ = OrderedDict([("mins", mins), ("maxs", maxs), ("nbinss", nbins), ("binarrays", binarrays), ("name", name), ("title", title)])

        self.TH1args, self.ndim = parrseHistArgs(mins,maxs,nbins,binarrays,name, title)
        self.__root_base_cls__=self.super = [ROOT.TH1F, ROOT.TH2F, ROOT.TH3F][self.ndim - 1]
        self.name = self.TH1args["name"]
        self.title = self.TH1args["title"]

        try:
            self.__root_base_cls__.__init__(self, *self.TH1args.values())
        except:
            print "Call to ROOT.TH{0} failed. Here are the arguments: {1}".format(self.ndim, self.TH1args.values())
            raise

        self.__ROOTAxis__ = [self.GetXaxis(),self.GetYaxis(),self.GetZaxis()]

        for _title, a in zip("xyz",self.__ROOTAxis__):
            a.SetTitle(_title)

        self.D3binvalues=np.zeros(tuple(map(lambda x: x.GetNbins(),self.__ROOTAxis__)))

        self.shape = tuple(map(lambda x:x.GetNbins(),self.__ROOTAxis__[0:self.ndim]))

        self.__axis_indicies__ = np.array(list(product(*tuple(map(lambda x: range(1,x+1),self.shape)))), np.int)

        self.binvalues = np.zeros(self.shape)
        self.binerrors = np.zeros(self.shape)


        self.bincenters=[[A.GetBinCenter(i) for i in range(1,A.GetNbins()+1)] for A in self.__ROOTAxis__[0:self.ndim]]
        self.__binLeftEdges__=tuple([np.array([A.GetBinLowEdge(i) for i in range(1,A.GetNbins()+1)]) for A in self.__ROOTAxis__[0:self.ndim]])
        self.__binRightEdges__=tuple([np.array([A.GetBinUpEdge(i) for i in range(1,A.GetNbins()+1)]) for A in self.__ROOTAxis__[0:self.ndim]])
        widths = [r - l for l, r in zip(self.__binLeftEdges__, self.__binRightEdges__)]

        def __product__(args):
            r = 1
            for n in args:
                r *= n
            return n
        # array of lengths, areas, or volumes of each bin with the same shape as self.binvalues.
        self.__bin_sizes__ = np.array(map(__product__, list(product(*widths)))).reshape(self.shape)

        self.binargs = [np.array([A.GetBinLowEdge(i) for i in range(1,A.GetNbins()+2)]) for A in self.__ROOTAxis__[0:self.ndim]]

        self.nbinsx = self.__ROOTAxis__[0].GetNbins()
        self.nbinsy = self.__ROOTAxis__[1].GetNbins()
        self.nbinsz = self.__ROOTAxis__[2].GetNbins()


        self.nentries=0

        self.nEntries = 0
        self.__needs_updating__ = False

        oubuilder = [range(1,len(i)+1) for i in self.bincenters]
        self.OverUnderFlows = [product(*(oubuilder[0:i]+oubuilder[i+1:])) for i in range(self.ndim)]

    def __copy__(self):
        result = [TH1F, TH2F, TH3F][self.ndim-1](**self.__self_args__)
        result.binvalues = self.binvalues.copy()
        result.binerrors = self.binerrors.copy()
        for i in range(3):
            result.__ROOTAxis__[i].SetTitle(self.__ROOTAxis__[i].GetTitle())
            result.SetTitleOffset(self.GetTitleOffset("XYZ"[i]))

        result.nEntries = self.nEntries
        result.SetEntries(self.nEntries)
        result.__update_hist_from_containers__()

        return result

    def __del__(self):
        pass

    def Fill(self, *args):
        self.nEntries+=1
        self.__needs_updating__=True
        self.__root_base_cls__.Fill(self,*args)

    def SetTitle(self,title):
        title = str(title)
        self.title = title
        self.super.SetTitle(self,title)

    def Draw(self, options = "", make_new_canvas = True,  *args):
        if "same" not in options.lower() and make_new_canvas:
            hist.tcanvii_refs.append(ROOT.TCanvas())
            hist.tcanvii_refs[-1].cd(0)
            hist.tcanvii_refs[-1].SetTitle("c"+str(len(hist.tcanvii_refs)-1))

        self.__root_base_cls__.Draw(self, options, *args)

    def __getArrayI_HistI__(self):
        return izip(product(*tuple(map(lambda x: xrange(1,x+1),self.shape))),product(*tuple(map(lambda x: xrange(0,x),self.shape))))

    def update_bin_containers_from_hist(self):

        for indicies4ROOT,indicies4Array in self.__getArrayI_HistI__():
            self.binvalues[indicies4Array] = self.GetBinContent(*indicies4ROOT)
            self.binerrors[indicies4Array] = self.GetBinError(*indicies4ROOT)

    def __update_hist_from_containers__(self):
        prev_entries = self.GetEntries()
        for indicies4ROOT,indicies4Array in self.__getArrayI_HistI__():
            bin_value = self.binvalues[indicies4Array]
            if np.isnan(bin_value) or np.isinf(bin_value):
                bin_value = 0
            args = indicies4ROOT+(bin_value,)
            self.super.SetBinContent(self,*args)

            bin_error = self.binerrors[indicies4Array]
            if np.isnan(bin_error) or np.isinf(bin_error):
                bin_error = 0
            args = indicies4ROOT + (bin_error,)
            self.super.SetBinError(self,*args)

        self.SetEntries(prev_entries)

    def __operator_decorator__(operator):

        def wrapper(self, *args, **kwargs):
            assert isinstance(self,hist)

            if self. __needs_updating__ or self.super.GetEntries(self)>self.nEntries:

                self.update_bin_containers_from_hist()
                self.__needs_updating__ = False

            return_value = operator(self,*args,**kwargs)

            good_bins = np.where(np.isnan(self.binvalues+self.binerrors)+np.isinf(self.binvalues + self.binerrors), False, True)

            self.binvalues = np.where(good_bins, self.binvalues, 0)
            self.binerrors = np.where(good_bins, self.binerrors, 0)

            self.__update_hist_from_containers__()

            return return_value

        return wrapper

    def __apply_method_wraper__(self,func):
        """
        :param func: A callable object that accepts ndim args, or an np.array and returns a scalar.
        :return: A callable object that accepts a tuple of lenth ndim, to be used in self.__apply__.
        """
        func_n_args = len(inspect.getargspec(func).args)
        if func_n_args == self.ndim:
            func_test = func(*((0.0,)*self.ndim))
            assert isinstance(func_test, numbers.Number), "Functions used in apply must return a scalar! not {0}".format(
                type(func_test))

            def wrapper(args):
                return func(*args)

        elif func_n_args == 1:
            func_test = func(np.array([0.0]*self.ndim))
            assert isinstance(func_test, numbers.Number), "Functions used in apply must return a scalar! not {0}".format(
                type(func_test))
            def wrapper(args):
                return func(np.array(list(args)))

        else:
            assert False, "When applying a function to a histogram, the function must ether take a tuple of length ndim, " \
                          "or a single np.array"
        return wrapper

    def __apply__(self, func):
        """ apply a function over all bin centers. returns result."""
        func  = self.__apply_method_wraper__(func)
        result = np.array(map(func, product(*self.bincenters))).reshape(self.shape)

        return result

    @__operator_decorator__
    def Add(self,other, c=1, error=None, princ_axis=None):

        if isinstance(other,numbers.Number):
            if error is not None:
                assert isinstance(error, numbers.Number), "The error must be a scalar in Add!"
                self.binerrors = np. sqrt(self.binerrors**2 + error**2)

            self.binvalues += c * np.ones(self.shape, dtype=np.float32) * other

        elif isinstance(other, hist):
            self.binerrors = np.sqrt(self.binerrors**2+(c*other.binerrors)**2)
            self.binvalues += c * other.binvalues

        elif princ_axis is None and isinstance(other, np.ndarray):
            assert self.shape == other.shape, "When not specifying 'princ_axis' in Add, the shape of array must match shape of histogram! Array: {0}".format(other)

            if error is not None:
                if isinstance(error,numbers.Number):
                    self.binerrors = np.sqrt(self.binerrors**2 + (error*np.ones(self.shape))**2)

                elif isinstance(error, np.ndarray):
                    assert self.shape == error.shape, "errors must of same shape!"
                    self.binerrors = np.sqrt(self.binerrors ** 2 + error ** 2)

                else:
                    raise NotImplementedError, "Error must be a scalar or an numpy array. "

            self.binvalues += other

        elif princ_axis is not None and isinstance(other, np.ndarray):
            if error is not None:
                if not isinstance(other, np.ndarray):
                    other = np.array(other)
                assert error.shape == other.shape, "errors are not of the correct shape!"
                self.binerrors = apply_operator_along_axis(self.binerrors, other,princ_axis, lambda a,b: np.sqrt(a**2+b**2))

            self.binvalues = apply_operator_along_axis(self.binvalues,other,princ_axis, operator=lambda a,b: a + b)


        elif hasattr(other,"__call__"):
            self.binvalues += self.__apply__(other)
            # TODO: Errors

    @__operator_decorator__
    def Multiply(self, other, error=None, princ_axis=None):
        if isinstance(other,numbers.Number):

            if error is not None:
                assert isinstance(error,numbers.Number), "<error> must be a scalar when multiplying by a scalar!"
                self.binerrors = np.sqrt((other*self.binerrors)**2 + (error*self.binvalues)**2)
            else:
                self.binerrors *= abs(other)

            self.binvalues *= other

        elif isinstance(other, hist):
            if self.shape == other.shape:
                self.binerrors = np.sqrt((other.binvalues*self.binerrors)**2 + (other.binerrors*self.binvalues)**2)
                self.binvalues *= other.binvalues

            elif  other.binvalues.shape != self.binvalues.shape:
                assert princ_axis is not None, "*** Must specify a principal axis when multiply {0}dim hist by {1}dim hist! *** {2}".format(
                    self.ndim, other.ndim, princ_axis)
                self.binerrors = np.sqrt(apply_operator_along_axis(self.binerrors, other.binvalues, princ_axis)**2 + apply_operator_along_axis(self.binvalues, other.binerrors, princ_axis)**2)
                self.binvalues = apply_operator_along_axis(self.binvalues, other.binvalues, princ_axis)

        elif isinstance(other, (np.ndarray, list, tuple)):
            other = np.array(other, dtype=np.float64)
            assert princ_axis is not None, "*** Must specify a principal axis when multiply {0}dim hist by {1}dim array! ***".format(self.ndim,other.ndim)
            if error is not None:
                assert other.shape == error.shape, "errors are of the wrong shape!"
                self.binerrors = np.sqrt(apply_operator_along_axis(self.binerrors,other,princ_axis)**2 + apply_operator_along_axis(self.binvalues,error,princ_axis)**2)
            else:
                self.binerrors = apply_operator_along_axis(self.binerrors, other, princ_axis = princ_axis, operator=lambda a,b:a*b)
            self.binvalues = apply_operator_along_axis(self.binvalues, other, princ_axis = princ_axis, operator=lambda a,b:a*b)

        elif hasattr(other, "__call__"): # TODO: have this function compute errors.

            if princ_axis is None:
                self.binvalues *= self.__apply__(other)
            else:
                assert len(princ_axis)==1, "2D arrays not supported for this operation "
                assert error is None, "The error argument not supported for this operation"
                assert isinstance(princ_axis, str)
                princ_axis = princ_axis.lower()
                assert princ_axis in "xyz", "invalid princ_axis"

                dummy_array = np.array(map(other, self.bincenters["xyz".index(princ_axis)]))
                self.Multiply(dummy_array, princ_axis = princ_axis)
        else:
            raise NotImplementedError

    @__operator_decorator__
    def Divide(self, other, const_err = 0):
        if isinstance(other, hist):

            old_settings = np.seterr(all='ignore')
            self.binerrors = np.sqrt((other.binerrors ** 2 * self.binvalues** 2 + self.binerrors ** 2 * other.binvalues** 2) / (other.binvalues**4))
            self.binvalues /= other.binvalues
            np.seterr(**old_settings)

        elif isinstance(other,numbers.Number):
            assert other!=0, "Divided by zero!"
            other = float(other)
            if const_err!=0:
                self.binerrors = np.sqrt((const_err ** 2 * self.binvalues ** 2 + self.binerrors ** 2 * other ** 2) / (
                        other ** 4))
            else:
                self.binerrors /= other

            self.binvalues /= other

        else:
             # TODO
            assert False, "TODO"

    def __ne__(self):
        copied_hist = self.__copy__()
        copied_hist.Multiply(-1)
        return copied_hist

    def __div__(self, other):
        copied_hist = self.__copy__()
        copied_hist.Divide(other)
        return copied_hist

    def __rdiv__(self, other):
        if isinstance(other, numbers.Number):
            other = float(other)
            result_hist = self.__copy__()
            result_hist.binerrors = np.zeros_like(self.binerrors)
            result_hist.binvalues = other * np.ones_like(self.binvalues)
            result_hist.__update_hist_from_containers__()
            result_hist.Divide(self)

        elif isinstance(other, hist):
            result_hist = other.__copy__()
            result_hist.Divide(self)

        elif isinstance(other, (np.ndarray,list)):
            other = np.ndarray(other, dtype=np.float32)
            assert self.binvalues.shape == other.shape, "rdiv called with incompatible shapes!"
            result_hist = self.__copy__()
            result_hist.binerrors = np.zeros_like(self.binerrors)
            result_hist.binvalues = other
            result_hist.__update_hist_from_containers__()
            result_hist.Divide(self)
        else:
            raise  NotImplementedError

        return result_hist

    @staticmethod
    def __getTGraph__(TH1,func, xerr_scale=1, y_errors=True, remove_zero_bins=False):
        """
        :param self:
        :param func: A function which takes three arguments, one corresponding to each axis, in the same order in which
         they appear.
        :param xerr_scale: Scale x errors if they are annoying.
        :param y_errors: Use y-errors?
        :param remove_zero_bins:
        :return: form: (TGraph, str)
        """
        TGraph_args, bin_info = TH1.__reduce_to_1D_by_scalar_function__(func, xerr_scale, y_errors, remove_zero_bins)
        graph = ROOT.TGraphErrors(*TGraph_args)

        graph.SetTitle(TH1.title)

        return graph, bin_info

    def normalize(self,width = True, squarIntegral = False):
        if not squarIntegral:
            if width:
                norm = abs(self)
            else:
                norm = np.sum(np.ndarray.flatten(self.binvalues))
        else:
            if width:
                norm = np.sum(np.ndarray.flatten(self.__bin_sizes__ * self.binvalues**2))
            else:
                norm = np.sum(np.ndarray.flatten(self.binvalues**2))

        if norm==0:
            warnings.warn("Trying to normalize histogram with zero integral!! Hist title: {0}".format(self.title))
            return 0
        else:
            self.Divide(norm)
            return norm

    def getErrorsHistogram(self, relative_error = True):
        newHist = self.__copy__()
        newHist.binvalues = self.binerrors
        if relative_error:
            newHist.binvalues /= self.binvalues
        newHist.binvalues = np.abs(newHist.binvalues)
        newHist.binerrors = np.zeros_like(self.binerrors)
        newHist.__update_hist_from_containers__()

        return newHist

    def transpose(self):
        assert self.ndim==2 and len(self.bincenters[0]) == len(self.bincenters[1])
        self.binvalues =  np.transpose(self.binvalues)
        self.binerrors = np.transpose(self.binerrors)
        self.__update_hist_from_containers__()
        return self

    def __idiv__(self, other):
        self.Divide(other)
        return self

    def __mul__(self, other):
        copied_hist = self.__copy__()
        copied_hist.Multiply(other)
        return copied_hist

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        self.Multiply(other)
        return self

    def __add__(self, other):
        copied_hist = self.__copy__()
        copied_hist.Add(other)
        return copied_hist

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        self.Add(other)
        return self

    def __sub__(self, other):
        copied_hist = self.__copy__()
        copied_hist.Add(other,-1)
        return copied_hist

    def __rsub__(self, other):
        copied_hist = self.__copy__()
        copied_hist = copied_hist.__ne__()  # h = -h
        copied_hist.Add(other)
        return copied_hist

    def __isub__(self, other):
        self.Add(other, -1)
        return self

    def __pow__(self, power, modulo=None):
        assert isinstance(power, numbers.Number) and power!=1
        copied_hist = self.__copy__()
        copied_hist.binvalues  = copied_hist.binvalues**power
        copied_hist.binerrors  = copied_hist.binerrors * copied_hist.binvalues * copied_hist.binvalues**(power-1)
        return copied_hist

    @__operator_decorator__
    def __abs__(self):
        widths_factor = 1
        if self.ndim >= 1:
            bx = self.__binRightEdges__[0] - self.__binLeftEdges__[0]
            if self.ndim == 1:
                widths_factor = bx

            elif self.ndim >= 2:
                by = self.__binRightEdges__[1] - self.__binLeftEdges__[1]
                if self.ndim == 2:
                    widths_factor = np.einsum('a,b->ab', bx, by)

                elif self.ndim == 3:
                    bz = self.__binRightEdges__[2] - self.__binLeftEdges__[2]
                    widths_factor = np.einsum('a,b,c->abc', bx, by, bz)

        norm = np.sum(np.ndarray.flatten(self.binvalues * widths_factor))
        return norm


    def reset(self):
        self.binvalues = np.zeros_like(self.binvalues)
        self.binerrors = np.zeros_like(self.binvalues)
        self.__update_hist_from_containers__()

    def __reduce_to_1D_by_scalar_function__(self,func, x_errors = False, y_errors = True, remove_zero_bins = True):

        results = self.__apply__(func)

        x = np.ndarray.flatten(results)
        y = np.ndarray.flatten(self.binvalues)
        erry = np.ndarray.flatten(self.binerrors)

        point_bin_correspondence = zip(zip(x,y,erry), product(*self.bincenters))
        bin_bininfo = ""

        srter = np.argsort(x)
        point_bin_correspondence = [point_bin_correspondence[_I_] for _I_ in srter]

        last_x = point_bin_correspondence[0][0][0]

        bin_bininfo = {"hist_bin_centers":[], "graph_points":[], "graph_errors":[]}

        for (_x_,_y_,_erry_), bin_center in point_bin_correspondence:
            if (round(last_x,2)!=round(_x_,2)):
                last_x = _x_

            bin_bininfo["hist_bin_centers"].append(bin_center)
            bin_bininfo["graph_errors"].append(_erry_)
            bin_bininfo["graph_points"].append((_x_,_y_))

        # TODO: fix errorx calculation
        errx = np.zeros_like(x)

        erry = np.ndarray.flatten(self.binerrors)


        if all(map(lambda x:x==0, y)):
            warnings.warn("All values in histogram are zero!")

        if  remove_zero_bins:
            selection = np.where(y!=0)
            x = x[selection]
            y = y[selection]
            errx = errx[selection]
            erry = erry[selection]

        if not y_errors:
            erry *= 0

        TGraph_args = (len(x),x,y,errx,erry)

        return TGraph_args, bin_bininfo

    @__operator_decorator__
    def set_max_rel_err(self, rel_err):
        old_settings = np.seterr(all='ignore')
        bin_rel_errs = self.binerrors/self.binvalues

        self.binvalues = np.where(bin_rel_errs < rel_err, self.binvalues,0 )

        self.binvalues = np.where(np.isnan(bin_rel_errs), 0,  self.binvalues)
        np.seterr(**old_settings)

        assert self.binvalues.shape==self.shape

    def Project(self, tree, _str_, cut = "", max_events = None, options= "", weight = None, start= None ):
        assert isinstance(_str_, str)

        __str__ = list(reversed(_str_.split(":")))

        assert len(list(__str__))==self.ndim or _str_=="" , "{0} does not match ndim = {1}".format(_str_, self.ndim)

        __str__ = ":".join(__str__)

        dummy_hist = self.__copy__()
        if max_events is not None:
            result = tree.Project(dummy_hist.name, __str__, cut, options, int(max_events), 0 if start is None else start)
        else:
            result = tree.Project(dummy_hist.name, __str__, cut, options)

        dummy_hist.update_bin_containers_from_hist()
        self.nEntries = self.GetEntries() + dummy_hist.GetEntries()
        self.SetEntries(self.nEntries)


        if weight is not None:
            assert isinstance(weight, numbers.Number)
            dummy_hist.Multiply(weight)

        self.Add(dummy_hist)

        return result


class TH1F(hist,ROOT.TH1F):

    def __init__(self,mins=None, maxs = None, nbinss = None, binwidths = None, binarrays = None, name=None, title=None):
        parsed_args = __parse_multidim_args_befor_hist__(1,mins,maxs,nbinss,binwidths,binarrays)+(name,)+(title,)
        hist.__init__(self,*parsed_args)

    def GetTGraph(self, x_errors = True):
        if x_errors:
            xerr = np.array(self.__binRightEdges__[0],dtype=np.float64)-np.array(self.__binLeftEdges__[0],dtype=np.float64)
        else:
            xerr = np.zeros_like(self.binvalues)
        graph = ROOT.TGraphErrors(len(self.binvalues),np.array(self.bincenters[0],dtype=np.float64),np.array(self.binvalues,dtype=np.float64),xerr, self.binerrors)
        graph.GetXaxis().SetTitle(self.GetXaxis().GetTitle())
        graph.GetYaxis().SetTitle(self.GetYaxis().GetTitle())
        return graph



class TH2F(hist,ROOT.TH2F):
    def __init__(self,mins=None, maxs = None, nbinss = None, binwidths = None, binarrays = None, name=None, title=None):
        parsed_args = __parse_multidim_args_befor_hist__(2,mins,maxs,nbinss,binwidths,binarrays)+(name,)+(title,)
        hist.__init__(self,*parsed_args)

    def get_1D_slices(self,fixed_axis, title_of_axis=None):
        assert isinstance(fixed_axis,str), "<fixed_axis> must be string!"
        fixed_axis = fixed_axis.lower()
        assert fixed_axis in ["x", "y" ], "ERROR: <fixed_axis> can be either 'x'or 'y', not {}".format(fixed_axis)

        if fixed_axis=="x":
            xy_index = 0
            array_of_bin_values = self.binvalues
            array_of_bin_errors = self.binerrors
            bin_args = self.binargs[1]

        elif fixed_axis=="y":
            xy_index = 1
            array_of_bin_values = np.transpose( self.binvalues)
            array_of_bin_errors = np.transpose( self.binerrors)
            bin_args = self.binargs[0]

        else: assert False  # should not be possible

        histos = []
        fixed_axis_values = []

        if title_of_axis is None:
            title_of_axis = ["x","y"][xy_index]

        for _index_ in range(len(array_of_bin_errors)):
            left_edge =  self.__binLeftEdges__[xy_index][_index_]
            right_edge =  self.__binRightEdges__[xy_index][_index_]
            hist = TH1F(binarrays =  bin_args, title= "{0} <= {1} < {2}".format(round(left_edge,2),title_of_axis,round(right_edge,2)))

            hist.binvalues = array_of_bin_values[_index_]
            hist.binerrors = array_of_bin_errors[_index_]
            hist.__update_hist_from_containers__()

            histos.append(hist)
            fixed_axis_values.append((left_edge,right_edge))

        return histos, fixed_axis_values

    def getTGraph(self,func, xerr_scale=1, y_errors=True, remove_zero_bins=False):
        return hist.__getTGraph__(self,func, xerr_scale, y_errors, remove_zero_bins)


class TH3F(hist,ROOT.TH3F):
    def __init__(self, mins=None, maxs=None, nbinss=None, binwidths=None, binarrays=None, name=None, title=None):
        parsed_args = __parse_multidim_args_befor_hist__(3, mins, maxs, nbinss, binwidths, binarrays) + (name,) + (
        title,)
        hist.__init__(self, *parsed_args)

    def get_2D_slices(self, fixed_axis, title_of_axis=None):
        assert isinstance(fixed_axis, str), "ERROR: <fixed_axis> must be a str!"
        fixed_axis = fixed_axis.lower()
        assert fixed_axis in ["x", "y", "z"], "ERROR: <fixed_axis> must be exactly 'x', 'y' , or 'z'!"
        tr_pose = [0, 1, 2]
        firsti = "xyz".index(fixed_axis)
        tr_pose[firsti], tr_pose[0] = tr_pose[0], tr_pose[firsti]
        tr_pose = tuple(tr_pose)

        array_of_bin_values = np.transpose(self.binvalues, tr_pose)
        array_of_bin_errors = np.transpose(self.binerrors, tr_pose)

        histos = []
        fixed_axis_values =[]

        if title_of_axis is None:
            title_of_axis = fixed_axis

        bin_args = [self.binargs[i] for i in tr_pose[1:]]

        for _index_ in range(len(array_of_bin_errors)):

            axis_right = self.__binRightEdges__[tr_pose[0]][_index_]
            axis_left = self.__binLeftEdges__[tr_pose[0]][_index_]
            new_TH2F = TH2F(binarrays=bin_args,title="{1} < {0} < {2}".format(title_of_axis, round(axis_left,2), round(axis_right,2)))
            new_TH2F.binvalues = array_of_bin_values[_index_]
            new_TH2F.binerrors = array_of_bin_errors[_index_]
            new_TH2F.__update_hist_from_containers__()
            histos.append(new_TH2F)
            fixed_axis_values.append((axis_left,axis_right))

        return histos, fixed_axis_values

    def getTGraph(self,func, xerr_scale=1, y_errors=True, remove_zero_bins=False):
        return hist.__getTGraph__(self,func, xerr_scale, y_errors, remove_zero_bins)


def get_from_ROOT_hist(hist):
    assert isinstance(hist,ROOT.TH1)
    dim = hist.GetDimension()

    result = None

    title = hist.GetTitle()
    if dim>=1:
        x_bins = map(hist.GetBinLowEdge, range(1,hist.GetXaxis().GetNbins()+2))
        if dim == 1:
            result = TH1F(binarrays=x_bins, title = title)
        if dim>=2:
            y_bins = map(hist.GetBinLowEdge, range(1,hist.GetYaxis().GetNbins()+2))
            if dim==2:
                result = TH2F(binarrays=[x_bins,y_bins], title=title)

            if dim>=3:
                z_bins = map(hist.GetBinLowEdge, range(1,hist.GetZaxis().GetNbins()+2))
                if dim == 3:
                    result = TH3F(binarrays=[x_bins,y_bins,z_bins], title=title)

    assert result is not None
    return result






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

