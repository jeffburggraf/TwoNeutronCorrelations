import re
import warnings
import numpy as np
import time
import itertools
import os
import sys

try:
    import ROOT
except:
    pass

def custom_formatwarning(msg, *a):
    # ignore everything except the message
    return "Warning: " + str(msg) + '\n'


warnings.formatwarning = custom_formatwarning

text_Padding_Width = 10

class branch_value(dict):

    def __init__(self, value=None, dtype='float', allow_zero_value=True, title="title", text_format="float"):

        self.title = title

        super(branch_value, self).__setitem__(0, value)
        self.allow_zero_value = allow_zero_value

        if dtype == "float":
            self.dtype = float
        elif dtype == "int":
            self.dtype = int
        else:
            assert False

        self.br = ROOT.vector(self.dtype)()

        if value is not None:
            self.br.push_back(self.dtype(value))
            self.__is_empty__ = False
        else:
            self.__is_empty__ = True

        # rest is for writing values to text file.
        self.integers_only = False
        self.none_text = "{{:<{0}d}}".format(text_Padding_Width).format(0)
        self.text_format = text_format

        self.__set_number_format__(text_Padding_Width)

        self.max_fitable_value = 10**(text_Padding_Width-1)

    def __set_number_format__(self, text_padding):
        if self.text_format == "float" and self.text_format != "int":
            self.number_format = "{{:<{0}.5f}}".format(text_padding)
        else:
            self.number_format = "{{:<{0}d}}".format(text_padding)
            self.integers_only = True

    def get_text(self):
        value = self[0]

        if value is None:
            return self.none_text

        if value > self.max_fitable_value:
            self.max_fitable_value*=10
            global text_Padding_Width
            text_Padding_Width +=1
            self.__set_number_format__(text_Padding_Width)


        if self.integers_only:
            value = int(value)

        return self.number_format.format(value)

    def __nonzero__(self):
        return self[0] is not None

    def __setitem__(self, key, value):
        assert key == 0

        if (value is None) or (value == 0 and not self.allow_zero_value):
            # self.__is_empty__ is used for optimization, since br.clear() is expensive.
            if not self.__is_empty__:  # if container is not empty, empty it.
                self.br.clear()
                self.__is_empty__ = True
                super(branch_value, self).__setitem__(key, None)
        else:
            super(branch_value, self).__setitem__(key, value)
            if self.br.size() == 0:
                self.br.push_back(np.float32(value))
            else:
                self.br[0] = value

            if self.__is_empty__:
                self.__is_empty__ = False

    def __getitem__(self, item):
        # print "getitem", id(self)
        return super(branch_value, self).__getitem__(item)

    def __repr__(self):
        return str(self[0])

    def __float__(self):
        if self:
            return float(self[0])
        else:
            return 0

    def __int__(self):
        if self:
            return int(self[0])
        else:
            return 0


class empty_contanor:
    def __init__(self):
        self.title = "empty"

    def __setitem__(self, key, value):
        pass

    def __getitem__(self, item):
        pass

    def __repr__(self):
        return "empty"


Empty = empty_contanor()


class event:
    """
    This class serves as the data container for each event.
    """

    def __init__(self):

        self.dE_warnings = 0

        self.next_event = branch_value(9000., title="next_event")
        self.nps = branch_value(title="nps", text_format="int")

        self.evt_type = branch_value(allow_zero_value=False, title="evt_type", text_format="int")
        self.cell = branch_value(allow_zero_value=False, title="cell", text_format="int")
        self.zaid = branch_value(allow_zero_value=False, title="zaid", text_format="int")
        self.mt = branch_value(allow_zero_value=False, title="mt", text_format="int")
        self.ter = branch_value(allow_zero_value=False, title="ter", text_format="int")
        self.sur = branch_value(allow_zero_value=False, title="sur", text_format="int")
        self.subsurf = branch_value(allow_zero_value=False, title="subsurf",  text_format="int") # Needed becasue of rounding error in surfaces such as 10.1, so 10.1 goes to sur=10 and subsurf = 1.
        self.par = branch_value(allow_zero_value=False, title="par", text_format="int")
        self.surf_theta = branch_value(allow_zero_value=False, title="surftheta",
                                       text_format="int")  # This value in the agle w.r.t. the surface in degrees (integer).
        self.bnk = branch_value(title="bnk", text_format="int")

        self.x = branch_value(title="x")
        self.y = branch_value(title="y")
        self.z = branch_value(title="z")
        self.dirx = branch_value(title="dirx")
        self.diry = branch_value(title="diry")
        self.dirz = branch_value(title="dirz")
        self.erg = branch_value(allow_zero_value=False, title="erg")
        self.wgt = branch_value(allow_zero_value=False, title="wgt")
        self.tme = branch_value(title="tme")
        self.SRCParID = branch_value(allow_zero_value=False,
                                     title="SRCParID")  # This value is incremented with each new source particle in a given NPS. Also is a TTree branch.

        self.dx = branch_value(allow_zero_value=False, title="dx")
        self.dE = branch_value(allow_zero_value=False, title="dE")

        self.previous_energy = None
        self.previous_position = [None, np.zeros((3,), np.float32)]

    def reset(self, soft=False):
        # If soft: reset everything except SRCParID, previous energy, and previous position.
        # Hard reset must be done each time the tracking of a new particle begins.

        if not soft:
            self.SRCParID[0] = None

            self.previous_energy = None
            self.previous_position[0] = None

        self.evt_type[0] = self.next_event[
                               0] // 1000  # Next event is subsequently updated in function read_next_event()

        self.dE[0] = None
        self.dx[0] = None

        self.cell[0] = None
        self.zaid[0] = None
        self.mt[0] = None
        self.ter[0] = None
        self.sur[0] = None
        self.subsurf[0] = None
        self.par[0] = None
        self.surf_theta[0] = None
        self.bnk[0] = None

        self.x[0] = None
        self.y[0] = None
        self.z[0] = None
        self.dirx[0] = None
        self.diry[0] = None
        self.dirz[0] = None
        self.erg[0] = None
        self.wgt[0] = None
        self.tme[0] = None

    def set_eE_dx(self, line):
        if self.evt_type[0] == 1. or self.evt_type[0] == 2.:
            self.previous_energy = self.erg[0]

            self.previous_position[1][0] = self.x[0]
            self.previous_position[1][1] = self.y[0]
            self.previous_position[1][2] = self.z[0]
            self.previous_position[0] = True
            self.dE[0] = None
            self.dx[0] = None

        elif self.evt_type[0] == 4.:
            if not self.previous_energy is None:
                if not self.erg[0] is None:
                    if self.mt[
                        0] != 18:  # In mcnp, neutrons seem to covertly change their identity during n_fiss events.
                        self.dE[0] = self.previous_energy - self.erg[0]
                    if self.dE[0] is not None and self.dE[0] < -0.0001:
                        if self.dE_warnings < 10:
                            warnings.warn(
                                "A particle increased its in energy between events (dE = {1})! See line {0} in txt file, if applicable".format(
                                    line, self.dE[0]))
                            self.dE_warnings += 1
                        elif self.dE_warnings == 10:
                            warnings.warn("Suppressing particle energy increase warnings.  ".format(line))
                            self.dE_warnings = 11

            if not self.previous_position[0] is None:
                self.dx[0] = np.linalg.norm(self.previous_position[1] - np.array([self.x[0], self.y[0], self.z[0]]))

            self.previous_position[1][0] = self.x[0]
            self.previous_position[1][1] = self.y[0]
            self.previous_position[1][2] = self.z[0]
            self.previous_energy = self.erg[0]

    def __Verbose_1__(self):
        event = str(int(self.evt_type))
        if self.bnk:
            event += "." + str(int(self.bnk[0]))
        else:
            event += "   "

        s = "nps: {0}  event: {1} mt: {2}".format(int(self.nps), event, int(self.mt) if self.mt else " ")
        return s

    def event_info(self, verbose=1):
        s = ""
        if verbose >= 1:
            s += self.__Verbose_1__()

        return s


Event = event()
PTRAC_BluePrints = {1.: [], 2.: [], 3.: [], 4.: [], 5.: []}


class event_writer():
    def __init__(self, path, write_level=2, POLIMI=False, write_text=False):
        assert isinstance(path, str)
        self.POLIMI = POLIMI

        self.line_number = 1

        if write_text:
            self.write_to_text_file = True
        else:
            self.write_to_text_file = False

        self.write_level = write_level
        assert write_level in [1, 2, 3], "write argument must be 1,2, or 3!"

        if path[-1] == "/":
            path = path[0:-1]

        assert os.path.isdir(path)
        self.path = path

        global ROOT

        try:
            import ROOT
            self.use_ROOT = True
        except:
            self.use_ROOT = False

        branches = []
        for _, attrib in Event.__dict__.iteritems():
            if isinstance(attrib, branch_value):
                branches.append(attrib)

        write_level_branch_to_remove = {1: ["surftheta", "sur","subsurf", "dx", "dE", "SRCParID", "bnk", "ter", "mt"],
                                        2: ["surftheta", "sur","subsurf", "dx", "dE", "SRCParID"],
                                        3: []}

        to_remove = ["next_event"]  # no need to make this a TTree branch.

        if not self.POLIMI:
            to_remove.append("SRCParID")

        for i in reversed(range(len(branches))):
            branch = branches[i]

            if branch.title in (to_remove + write_level_branch_to_remove[self.write_level]):
                del branches[i]

        if self.write_level != 3:
            msg = "Some available information is not being written to TTree.\n" \
                  "Values not being written: {}\n" \
                  "To write all values, call main() with keyword argument 'write=3'." \
                .format(write_level_branch_to_remove[write_level])

            warnings.warn(msg)

        ordered_branch_values = ["nps", "evt_type", "par", "bnk", "mt", "ter", "cell", "zaid", "erg", "x", "y", "z",
                                 "dirx", "diry", "dirz", "tme", "wgt", "dE", "dx", "sur","subsurf"]

        def sort_func(value):
            assert isinstance(value, branch_value)
            if value.title in ordered_branch_values:
                return ordered_branch_values.index(value.title)
            else:
                return len(ordered_branch_values)

        self.branches = sorted(branches, key=sort_func)

        # print "\n".join([b.title for b in self.branches])

        if not self.use_ROOT:
            warnings.warn("No ROOT module available. Data will not be placed into a ROOT TTree.")
            self.write_to_text_file = True

        if self.write_to_text_file:
            self.txt_file = open("{0}/parsed_ptrac.txt".format(self.path), "w")
            print "Creating text file at {0}/parsed_ptrac.txt\n".format(self.path)

            self.write_to_text_file = True

            _columb_label_ = "{{:{}s}}".format(text_Padding_Width)
            first_line = "".join([_columb_label_.format(b.title) for b in self.branches])
            self.txt_file.write(first_line + "\n")

        if self.use_ROOT:
            self.file = ROOT.TFile("{0}/ptrac.root".format(path), "recreate")
            self.tree = ROOT.TTree("tree", "tree")

            for b in self.branches:
                assert isinstance(b, branch_value)
                self.tree.Branch(b.title, b.br)

    def write_event(self):
        if Event.tme:
            Event.tme[0] *= 0.1  # shakes to ns.

        if Event.sur:
            sur = Event.sur[0]
            if int(sur)!=sur:
                Event.subsurf[0] = int(round((sur-int(sur))/0.1,1))

            Event.sur[0] = int(sur)

        if self.use_ROOT:
            self.tree.Fill()

        if self.write_to_text_file:
            self.line_number += 1
            self.txt_file.write("".join([b.get_text() for b in self.branches]))

            self.txt_file.write("\n")

    def finish(self):
        if self.use_ROOT:
            self.tree.Write()
            self.tb = ROOT.TBrowser()

        if self.write_to_text_file:
            self.txt_file.close()


def header(_file):
    # assert isinstance(cls, __reader__)

    assert isinstance(_file, file)

    path = _file.name

    assert re.match(" +-1", _file.readline()), "First line must be '     -1'.  Something is wrong. "

    _line_ = _file.readline()

    version_match = re.match("mcnp +([0-9]) +(.+)", _line_)

    if not version_match:
        version_match = re.match("mcnp([a-z])(.+)", _line_)

    version, date_time = None, None

    if version_match:
        version, date_time = version_match.groups()

    if version not in ["6", "5", "x"]:
        warnings.warn("This code has not been tested with this version of MCNP.")

    INP_title = re.sub(" ", "", _file.readline())

    print "Processing file '{0}' from input deck titled '{1}'\n".format(os.path.basename(path), INP_title[:-1])

    current_line = _file.readline()

    entry_iter = []

    last_pos = None

    _i_ = 0

    while True:
        if _i_ > 3:
            warnings.warn("There seems to be more than three variable ID lines, which is unusual. Be warned. ")

        _m = re.match(r" +\d\.\d+E\+\d+ ", current_line)

        if not _m:
            break
        else:
            entry_iter += [float(i) for i in current_line.split()]

        last_pos = _file.tell()
        current_line = _file.readline()
        _i_ += 1

    _file.seek(last_pos)
    i = 1
    keys = {}
    keyWordStrs = ["BUFFER", "CELL", "EVENT", "FILE", "FILTER", "MAX", "MENP", "NPS", "SURFACE", "TALLY", "TYPE",
                   "VALUE", "WRITE", "", "", ""]  # from MCNP6 user manual appendix F

    while i < len(entry_iter):
        num_key_entries = int(entry_iter[i])
        key_values = []
        if int(num_key_entries) != 0:
            for Di in range(1, int(num_key_entries) + 1):
                i = i + 1

                if not (i) < len(entry_iter):
                    break
                key_values.append(entry_iter[i])

        keys[keyWordStrs[len(keys)]] = key_values
        i += 1

    if keys["WRITE"][0] != 2.0:
        warnings.warn(
            "If the result seems incorrect, you might want to change the write keyword to 'all' on the PTRAC card (recomended!)")

    n_vars_line = map(lambda x: int(float(x)), _file.readline().split())

    if len(n_vars_line) != 20:
        warnings.warn(
            "This line should have 20 entries (see code), as per appendix F ('pg' 13-2) in the MCNP6 usermanual."
            "Not a big cause for concern...")

    n_var = [[n_vars_line[0]], n_vars_line[1:3], n_vars_line[3:5], n_vars_line[5:7], n_vars_line[7:9],
             n_vars_line[9:11]]

    if n_var[0] != [2]:
        warnings.warn("Setting the number of variables on an NPS line to 2, even though it is apparently {0}. "
                      "This value should always be 2!".format(n_var[0][0]))

    single_particle = n_vars_line[11]

    var_ID_lines = []

    num_M_lines = 2 if keys["WRITE"][0] != 2.0 else 3

    for i in range(num_M_lines):
        var_ID_lines += map(lambda x: int(float(x)), _file.readline().split())

    NPS_IDS = var_ID_lines[:n_var[0][0]]

    assert NPS_IDS == [1, 2], \
        "This algorithm assumes that the variable ids for the two entries on NPS line are always 1, and 2. "

    del var_ID_lines[:n_var[0][0]]

    SRC_IDS0 = var_ID_lines[:n_var[1][0]]
    del var_ID_lines[:n_var[1][0]]

    SRC_IDS1 = var_ID_lines[:n_var[1][1]]
    del var_ID_lines[:n_var[1][1]]

    BNK_IDS0 = var_ID_lines[:n_var[2][0]]
    del var_ID_lines[:n_var[2][0]]

    BNK_IDS1 = var_ID_lines[:n_var[2][1]]
    del var_ID_lines[:n_var[2][1]]

    SUR_IDS0 = var_ID_lines[:n_var[3][0]]
    del var_ID_lines[:n_var[3][0]]

    SUR_IDS1 = var_ID_lines[:n_var[3][1]]
    del var_ID_lines[:n_var[3][1]]

    COL_IDS0 = var_ID_lines[:n_var[4][0]]
    del var_ID_lines[:n_var[4][0]]

    COL_IDS1 = var_ID_lines[:n_var[4][1]]
    del var_ID_lines[:n_var[4][1]]

    TER_IDS0 = var_ID_lines[:n_var[5][0]]
    del var_ID_lines[:n_var[5][0]]

    TER_IDS1 = var_ID_lines[:n_var[5][1]]
    del var_ID_lines[:n_var[5][1]]

    # Print the lines below to see what the code determined were the variable ID's, in order, for each type of line.
    # print "\nNPS:  ",NPS_IDS
    # print "SRC0: ",SRC_IDS0
    # print "BNK0: ",BNK_IDS0
    # print "SUR0: ",SUR_IDS0
    # print "COL0: ",COL_IDS0
    # print "TER0: ",TER_IDS0

    # print "SRC1: ",SRC_IDS1
    # print "BNK1: ",BNK_IDS1
    # print "SUR1: ",SUR_IDS1
    # print "COL1: ",COL_IDS1
    # print "TER1: ",TER_IDS1

    if len(var_ID_lines) != 0:
        warnings.warn("PTRAC 'Variable ID lines' did not add up as expected!"" The same number of variables declared "
                      "for SRC,BNK, ect., should appear in the subsequent variable ID line. ")

    assert 7 == SRC_IDS0[0] == BNK_IDS0[0] == SUR_IDS0[0] == COL_IDS0[0] == TER_IDS0[0], \
        "This algorithm always assumes the 'next event' indicator " \
        "is always at the beginning of each event line, so these must all be equal to 7!"

    last_pos = _file.tell()

    global G_ptrac_begin_pos
    G_ptrac_begin_pos = last_pos

    first_event_line = _file.readline().split()

    if not first_event_line:
        assert False, "PTRAC file has no events!!"
    assert re.match("[0-9]+", first_event_line[0]).group(0) == first_event_line[0], \
        "The first entry should be an integer, namely, it's the first NPS of the PTRAC file. " \
        "If this check fails, something is wrong. "

    assert len(first_event_line) == 2, "Something is wrong. First event line is not of correct form. "
    assert len(first_event_line[1]) == 4, "Something is wrong. First event line is not of correct form. "

    _file.seek(G_ptrac_begin_pos)

    """
    The two dictionaries, index_dict1 and index_dict0 (below), are used here to link 'pointers' in python
    to the values as they appear in the PTRAC file. These dictionaries are used to construct a list
    of references to class attributes, which is sorted such that the order of the attributes matches the order of the
    variables are they appear on the two lines corresponding every event.
    Values of the numbers are taken from the MCNP user manual appendix F, table 13-4.  
    """

    index_dict1 = {20: "x", 21: "y", 22: "z", 23: "dirx", 24: "diry", 25: "dirz", 26: "erg", 27: "wgt", 28: "tme"}
    index_dict0 = {13: "surf_theta", 7: "next_event", 3: "cell", 17: "cell", 10: "zaid", 11: "mt",
                   12: "sur", 14: "ter", 16: "par"}  # TODO: why do 3 an 17 refer to cell variable?

    global PTRAC_BluePrints
    try:
        PTRAC_BluePrints[1.].append([(getattr(Event, index_dict0[index]) if index in index_dict0 else Empty) for
                                     index in SRC_IDS0])
        PTRAC_BluePrints[2.].append([(getattr(Event, index_dict0[index]) if index in index_dict0 else Empty) for
                                     index in BNK_IDS0])
        PTRAC_BluePrints[3.].append([(getattr(Event, index_dict0[index]) if index in index_dict0 else Empty) for
                                     index in SUR_IDS0])
        PTRAC_BluePrints[4.].append([(getattr(Event, index_dict0[index]) if index in index_dict0 else Empty) for
                                     index in COL_IDS0])
        PTRAC_BluePrints[5.].append([(getattr(Event, index_dict0[index]) if index in index_dict0 else Empty) for
                                     index in TER_IDS0])

        PTRAC_BluePrints[1.].append([(getattr(Event, index_dict1[index]) if index in index_dict1 else Empty) for
                                     index in SRC_IDS1])
        PTRAC_BluePrints[2.].append([(getattr(Event, index_dict1[index]) if index in index_dict1 else Empty) for
                                     index in BNK_IDS1])
        PTRAC_BluePrints[3.].append([(getattr(Event, index_dict1[index]) if index in index_dict1 else Empty) for
                                     index in SUR_IDS1])
        PTRAC_BluePrints[4.].append([(getattr(Event, index_dict1[index]) if index in index_dict1 else Empty) for
                                     index in COL_IDS1])
        PTRAC_BluePrints[5.].append([(getattr(Event, index_dict1[index]) if index in index_dict1 else Empty) for
                                     index in TER_IDS1])

        PTRAC_BluePrints[9.] = [[Event.nps, Event.next_event]]

    except AttributeError:
        assert False, "Attribute error! Every value in index_dict0 and index_dict1 must also be an attribute of class object 'event'."


class progress_mon:
    def __init__(self, _file_ ):
        assert isinstance(_file_, file)

        self.file = _file_
        self.size = float(os.stat(_file_.name).st_size)

        self.next_report = time.time() + 2
        self.bytes = []
        self.times = []
        self.n_events = []

    def progress(self,events):
        if  time.time() >= self.next_report:
            percent_complete = round(100* self.file.tell()/self.size,2)
            self.next_report =  time.time() + 2

            msg = "{0}% ".format(percent_complete)

            if len(self.bytes)!=3:
                self.bytes.append(self.file.tell())
                self.times.append(time.time())
                self.n_events.append(events)

            if len(self.bytes)==3:
                self.bytes[0:-1] = self.bytes[1:]
                self.bytes[-1] = self.file.tell()

                self.times[0:-1] = self.times[1:]
                self.times[-1] = time.time()

                self.n_events[0:-1] = self.n_events[1:]
                self.n_events[-1] = events

                bytes_per_second = (self.bytes[-1] - self.bytes[0])/(self.times[-1] - self.times[0])
                sec_remaining = int((self.size - self.file.tell())/bytes_per_second)
                hours = sec_remaining//(60*60)
                minutes = (sec_remaining%(60*60))//60
                sec = sec_remaining%(60)
                remaining = ""
                if hours:
                    remaining += "{0} hr ".format(hours)
                if minutes:
                    remaining += "{0} min ".format(minutes)
                if sec:
                    remaining += "{0} sec ".format(sec)

                if remaining:
                    msg += ", {0}remaining".format(remaining)

                event_rate = int((self.n_events[-1] - self.n_events[0]) / (self.times[-1] - self.times[0]))

                msg += ", {0} events/s".format(event_rate)


            print msg




# @profile

def read_next_event(evt,_file):
    """
    Reads and the next 2 lines, or only a single line if the file is positioned at an nps line, which occur once at the beginning of every
    particle history.
    After reading the line(s) in to an array(s) of floats, each float is placed into its respective container.
    The containers are an attribute of the global event object named Event. Each container is an instance of the
    branch_value class -- which is a user-defined python object used to automatically create and update a pointer (to a c++ std::vector)
    for a ROOT TTree leaf. The branch_value class uses the __setitem__ method to update its pointer, eg "branch_value[0]=value".

    :param evt: Global event obj, named Event
    :param file: open ptrac file
    :return: True if event processed, None if end of file is reached.
    """

    read_next_event.position_at_last_call = _file.tell()

    assert isinstance(evt, event) # not needed

    next_event = Event.next_event[0]

    try:
        blue_print = PTRAC_BluePrints[next_event//1000]
    except KeyError:
        print "\nError! Expected 'next_event' is not a valid event!\nnext_event: {0}\nvalid event types: {1}\nnext line to read:\n{2}"\
            .format(evt.next_event,[1000*i for i in PTRAC_BluePrints.keys()],_file.readline())
        time.sleep(0.5)
        assert 0

    evt.reset(soft=True)

    if Event.next_event[0]//1000==2.:
        Event.bnk[0] = Event.next_event[0] % 2000

    for refs in blue_print:
        line_floats = np.fromstring(_file.readline(), float, sep=" ")

        if len(refs) != len(line_floats):
            if len(line_floats)==0:
                return None
            else:
                print "line 'blueprint' mismatch! Unexpected number of elements on this line.\n " \
                      "Expected values on line: {0}\nline values read: {1}" \
                    .format([ref.title for ref in refs], line_floats)
                time.sleep(0.5)
                assert 0

        [ref.__setitem__(0, value) for ref, value in itertools.izip(refs, line_floats)]

    return True


POLIMI_src_msg_1= "A POLIMI source particle has appeared at more than one location! w/o a constant location for \n" \
                  "POLIMI source particles, banked photonuclear neutrons cannot be distinguished from POLIMI source neutrons,\n"\
                  "since both are, unfortunately, assigned to the same bank number."
Path = None

class __event_interpreter__():
    def __init__(self, write_level, using_POLIMI_source = False, write_untracked_POLIMI_src_events = False, write_to_text_file = False):
        self.POLIMI_SOURCE = using_POLIMI_source

        self.write_untracked_POLIMI_src_events = write_untracked_POLIMI_src_events

        self.write_level = write_level

        self.writer = event_writer(Path, write_level=write_level,POLIMI = using_POLIMI_source, write_text=write_to_text_file)

        self.SRC_events_seen_this_nps = 0

        self.error_check = True # Extra error checking is performed. Turned off after 1000 events are successfully read.

        self.events_processed = 0

        if self.POLIMI_SOURCE:
            self.polimi_tracking_begin = False
            self.last_src_pos = None
            self.POLIMI_source_pos = None

            # When using A POLIMI source, there is always a 14 MeV neutron at the beginning of each NPS that is
            # immediately terminated
            self.POLIMIfourteenMeVDefault = None

    def process_next(self):
        if not read_next_event(Event, _file):
            return None
        else:
            # check if an NPS line was just read.
            if Event.evt_type[0] == 9.:
                self.SRC_events_seen_this_nps = 0
                if self.POLIMI_SOURCE:
                    # reset POLIMI specific variables
                    self.polimi_tracking_begin = False

                Event.reset(soft=False) # End of nps, so values must be reset.

                # return false so that this event (NPS event) is not written.
                return False

            if self.POLIMI_SOURCE:
                # Special handling is needed until MCNP begins the actual tracking of POLIMI source particles (see note 1)TODO: note 1
                if not self.polimi_tracking_begin:
                    assert Event.evt_type[0] == 1.0, "Only SRC events should appear before POLIMI tracking begins." \
                                                       "Make sure your PTRAC card has SRC events written to the file!"

                    if self.write_untracked_POLIMI_src_events:
                        if Event.erg[0]!=14.:
                            warnings.warn("When using a POLIMI source, a default 14 MeV neutron event is expected to appear at the beginning of each NPS.\n")
                        self.polimi_tracking_begin = True

                    elif self.error_check<1000:
                        # Until error checking is done, check for the 14 MeV default neutron event at the beginning of each NPS.
                        if not self.POLIMIfourteenMeVDefault:
                            if Event.erg[0]!=14.:
                                warnings.warn("When using a POLIMI source, a default 14 MeV neutron event is expected to appear at the beginning of each NPS.\n")

                            self.POLIMIfourteenMeVDefault = True

                    # check if event is NOT a SRC event.
                    if Event.next_event[0] != 1000.:
                        # POLIMI has finished its sequence of bullshit SRC particles, so now real particle tracking begins.
                        self.polimi_tracking_begin = True

                        _file.seek(read_next_event.position_at_last_call) # seek file read position to last read event.

                        Event.SRCParID[0] = 1

                        Event.next_event[0] = 1000. # Because we just went back 1 event, the next event is a SRC event.

                        if self.error_check < 1000:
                            if self.POLIMI_source_pos is None:
                                self.POLIMI_source_pos = Event.x[0], Event.y[0], Event.z[0]
                            else:
                                if self.POLIMI_source_pos != (Event.x[0], Event.y[0], Event.z[0]):
                                    global POLIMI_src_msg_1
                                    if POLIMI_src_msg_1:
                                        warnings.warn(POLIMI_src_msg_1)
                                        POLIMI_src_msg_1 = None

                    return False  # return false because the tracking of POLIMI source particle hasn't begun yet.
                else:
                    if Event.evt_type[0] == 2.:
                        if Event.bnk[0] == 25 or Event.bnk[0] == 26:
                            if (Event.x[0], Event.y[0], Event.z[0]) == self.POLIMI_source_pos:
                                Event.bnk[0] = None
                                Event.evt_type[0] = 1.
                                Event.SRCParID[0] += 1

            if self.write_level>1:
                Event.set_eE_dx(self.writer.line_number)

            if self.error_check < 1000:
                self.error_check += 1

                if Event.evt_type[0] == 1.:
                    self.SRC_events_seen_this_nps += 1

                    if self.SRC_events_seen_this_nps > 1 and not self.POLIMI_SOURCE:
                        assert 0, "More than one SRC event seen in a single NPS. Turn POLIMI on! (This can also happen due to lost particles... Check for errors.)"

            self.events_processed += 1
            return True

def main(path,max_events = None, using_POLIMI_source=False, write_untracked_POLIMI_src_events = False, verbose = 2, write = 3, write_to_text_file = False):

    global _file,Path
    _file = open(path)

    Path = os.path.dirname(path)

    header(_file)

    global evt_interp

    evt_interp = __event_interpreter__(using_POLIMI_source=using_POLIMI_source,write_level = write, write_untracked_POLIMI_src_events = write_untracked_POLIMI_src_events, write_to_text_file = write_to_text_file)

    t1 = time.time()

    progress = progress_mon(_file)

    n_events = 0

    try:
        while True:
            has_written_current_event = False  # to write final event in case of error. Aids debugging.

            read_event_result = evt_interp.process_next()

            if read_event_result is None:
                break

            if read_event_result:
                evt_interp.writer.write_event()
                has_written_current_event = True
                n_events+=1

                progress.progress(n_events)

                if verbose:
                    Event.event_info(verbose)

                if max_events:
                    if n_events>max_events:
                        print "Terminating after max events reached..."
                        break
    except:
        # It is very useful to write the text file (if applicable) in the event of an error.
        if not has_written_current_event:
            evt_interp.writer.write_event()
        evt_interp.writer.finish()
        raise

    evt_per_sec = n_events/(time.time() - t1)
    #TODO check weight histogram!

    print "Mean events per second: {}".format(int(evt_per_sec))

    evt_interp.writer.finish()

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('path_to_ptrac',type=str)
    parser.add_argument('-max_events',type=int, default = None)
    parser.add_argument('-write',type=int,default=3)

    parser.add_argument('-use_POLIMI_src',            action='store_const', const=True, default=False)
    parser.add_argument('-write_all_POLIMI_src_evts', action='store_const', const=True, default=False)
    parser.add_argument('-write_text_file',           action='store_const', const=True, default=False)

    print parser.print_help()

