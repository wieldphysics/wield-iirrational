# $Id: foton.py 7372 2015-05-13 20:57:23Z james.batch@LIGO.ORG $
# Python module for foton scripting, allows generation of filter coefficients
# from foton design strings, and parsing of foton filter files.
# Author: Christopher Wipf, LIGO - Massachusetts Institute of Technology
# May 6, 2014


import sys
import os
import numpy as np
import warnings

try:
    import ROOT
except ImportError:
    from . import ROOT
import cppyy
import ctypes

from array import array

from .. import TFmath
from . import conversion

ROOT.gInterpreter.AddIncludePath("/usr/include/gds")
ROOT.gSystem.AddDynamicPath("/usr/lib64")
ROOT.gSystem.Load("libRdmtsigp")
ROOT.gSystem.Load("libRgdsplot")

# some tests of ROOT to cause an early failure at IMPORT if missing
ROOT.filterwiz.FilterFile
# ROOT.iirorder
ROOT.FilterDesign
ROOT.iir2z
ROOT.iir2zpk
ROOT.iir2direct
ROOT.iirsoscount
ROOT.iir2poly

TcplxD = ROOT.basicplx(ROOT.Double)
TcplxF = ROOT.basicplx(float)

__docformat__ = "restructuredtext"
__all__ = [
    "FilterFile",
    "iir2zpk",
    "iir2zpkstr",
]


def _reverse_dict(d):
    return dict((idx, val) for val, idx in d)


input_switches = {
    "AlwaysOn": 1,
    "always_on": 1,
    "alwayson": 1,
    "ZeroHistory": 2,
    "zero_history": 2,
    "zerohistory": 2,
}

input_switches_R = [
    "always_on",
    "zero_history",
]

output_switches = {
    "immediately": 1,
    "Immediately": 1,
    "immediate": 1,
    "Immediate": 1,
    "ramp": 2,
    "Ramp": 2,
    "inputcrossing": 3,
    "InputCrossing": 3,
    "input_crossing": 3,
    "zerocrossing": 4,
    "ZeroCrossing": 4,
    "zero_crossing": 4,
}
output_switches_R = [
    "immediate",
    "ramp",
    "input_crossing",
    "zero_crossing",
]

design_types = [
    "LowPass",
    "HighPass",
    "BandPass",
    "BandStop",
]

if False:
    # to calm lint on py2
    unicode = None

if sys.version_info > (3, 0):
    str_types = (str,)
else:
    str_types = (str, unicode)


def check_ftype(ftype, F1_Hz, F2_Hz):
    assert ftype in design_types
    assert F1_Hz is not None
    if ftype in ["BandPass", "BandStop"]:
        assert F2_Hz is not None
        assert F2_Hz > F1_Hz
        return True
    return False


class FilterFile(object):
    """
    Read and edit CDS filter files.

    Example:
      >>> from ligocds.foton import FilterFile
      >>> filterfile = FilterFile('/cvs/cds/mit/chans/M1PDE.txt')
      >>> filterfile['LSC_CARM'][0].rate
      65536.0
      >>> filterfile['LSC_CARM']['PHlead'].name
      'PHlead'

    """

    def __init__(
        self,
        filename=None,
    ):
        self._filter_file_raw = ROOT.filterwiz.FilterFile()
        self.filename = filename
        if filename:
            self.read(filename)

    def __getitem__(self, key):
        lookup = lambda: self._filter_file_raw.find(key)
        # The subtlety here is that find() returns a C pointer that
        # can be invalidated when modules are added or removed -- so a
        # 'Module' needs to carry around this lookup function to
        # retrieve a valid pointer.
        try:
            item = lookup()
        except:
            raise KeyError(key)
        # need to use specifically this syntax for ROOT to check Null pointers
        # The python standard is None does NOT WORK
        if item == None:
            raise KeyError(key)
        return Module(lookup)

    def __setitem__(self, key, val):
        if not isinstance(val, Module):
            raise ValueError(val)
        self._filter_file_raw.add(key, val.rate)
        for sec in val:
            self[key][sec.index] = sec

    def __delitem__(self, key):
        self._filter_file_raw.remove(key)

    def __contains__(self, key):
        try:
            item = self._filter_file_raw.find(key)
        except:
            return False
        if item is None:
            return False
        else:
            return True

    def __iter__(self):
        for key in self.keys():
            yield key

    def keys(self):
        return [fm.getName() for fm in self._filter_file_raw.modules()]

    def items(self):
        return [(key, self[key]) for key in self.keys()]

    def refresh(self):
        return self._filter_file_raw.update()

    def valid(self):
        val = True
        for name, fm in self.items():
            for sec in fm:
                if not sec.valid():
                    val = False
                    print(name + "[" + str(sec.index) + "]", file=sys.stderr)
                    print("invalid", file=sys.stderr)
        return val

    def read(self, filename):
        """
        Load foton filter file. If filename does not end with an extension and has no path separators,
        then it is assumed to be a system name, and it is loaded from
        /opt/rtcds/$SITE/$IFO/chans/filename.upper().txt
        """
        fbase, ext = os.path.splitext(filename)
        if not ext:
            a, b = os.path.split(filename)
            if not a:
                site = os.environ["site"]
                ifo = os.environ["ifo"]
                filename = "/opt/rtcds/{0}/{1}/chans/{2}.txt".format(
                    site, ifo, filename.upper()
                )
                print("Appears to be a system name, reading from")
                print(filename)
        self.filename = os.path.abspath(filename)
        self._filter_file_raw.read(self.filename)

    def write(self, *args):
        "Save filter file."
        if not self.filename:
            raise Exception("undefined filename")
        if not (self.valid and self.refresh() and self.valid):
            raise Exception("Filters Not All Valid!")
        if len(args) == 0:
            self._filter_file_raw.write(self.filename)
        else:
            self._filter_file_raw.write(*args)


class Module(object):
    def __init__(self, lookup):
        self._module_raw = lookup

    @property
    def fm(self):
        return self._module_raw()

    @property
    def name(self):
        return self.fm.getName()

    @name.setter
    def name(self, val):
        return self.fm.setName(val)

    @property
    def rate(self):
        return self.fm.getFSample()

    @rate.setter
    def rate(self, val):
        return self.fm.setFSample(val)

    def __len__(self):
        return ROOT.filterwiz.kMaxFilterSections

    def __getitem__(self, key):
        if isinstance(key, str_types):
            for n in range(len(self)):
                item = self.fm[n]
                if item.getName() == key:
                    return Section(self._module_raw, n)
            raise KeyError(key)
        elif key in range(len(self)):
            return Section(self._module_raw, key)
        else:
            raise KeyError(key)

    def __setitem__(self, key, val):
        self[key].copyfrom(val)
        if isinstance(key, str_types):
            self[key].name = key

    def __contains__(self, key):
        try:
            # access just to test for presence
            self[key]
            return True
        except:
            return False

    def __iter__(self):
        for n in range(len(self)):
            yield self[n]


class Section(object):
    def __init__(self, lookup_fm, key):
        self._lookup_fm = lookup_fm
        self._key = key
        if not self.valid:
            warnings.warn(
                "This module '{0}' does not currently have valid design string!".format(
                    self.name
                )
            )

    def xfer(self, F_Hz):
        return iir2xfer(self.filter_raw, F_Hz)

    @property
    def ZPKs(self):
        return iir2zpk(self.filter_raw, plane="s")

    @property
    def ZPKsf(self):
        return iir2zpk(self.filter_raw, plane="f")

    @property
    def ZPKf(self):
        return iir2zpk(self.filter_raw, plane="f")

    @property
    def ZPKz(self):
        return iir2zroots(self.filter_raw)

    @property
    def zroots(self):
        return iir2zroots(self.filter_raw)

    @property
    def SOS(self):
        return iir2sos(self.filter_raw)

    @property
    def sec(self):
        return self._lookup_fm()[self._key]

    @property
    def F_sample_Hz(self):
        return self.filter_raw.getFSample()

    @property
    def F_nyquist_Hz(self):
        return self.filter_raw.getFSample() / 2

    @property
    def index(self):
        return self.sec.getIndex()

    @index.setter
    def index(self, val):
        return self.sec.setIndex(val)

    @property
    def name(self):
        return self.sec.getName()

    @name.setter
    def name(self, val):
        return self.sec.setName(val)

    @property
    def design(self):
        return self.sec.getDesign()

    @design.setter
    def design(self, val):
        return self._set_design(val)

    @property
    def filter_raw(self):
        return self.sec.filter()

    @property
    def order(self):
        return ROOT.iirorder(self.filt.get())

    @property
    def input_switch(self):
        return input_switches_R[self.sec.getInputSwitch()]

    @input_switch.setter
    def input_switch(self, val):
        return self.sec.setInputSwitch(input_switches[val])

    @property
    def output_switch(self):
        return output_switches_R[self.sec.getOutputSwitch()]

    @output_switch.setter
    def output_switch(self, val):
        return self.sec.setOutputSwitch(output_switches[val])

    @property
    def ramp(self):
        return self.sec.getRamp()

    @ramp.setter
    def ramp(self, val):
        return self.sec.setRamp(val)

    @property
    def tolerance(self):
        return self.sec.getTolerance()

    @tolerance.setter
    def tolerance(self, val):
        return self.sec.setTolerance(val)

    @property
    def timeout(self):
        return self.sec.getTimeout()

    @timeout.setter
    def timeout(self, val):
        return self.sec.setTimeout(val)

    # @property
    # def header(self):
    #        return self.sec.getHeader()

    # @header.setter
    # def header(self, val):
    #        return self.sec.setHeader(val)

    @property
    def valid(self):
        return self.sec.valid()

    def empty(self):
        return self.sec.empty()

    def check(self):
        return self.sec.check()

    def refresh(self):
        return self.sec.update()

    def update(self):
        return self.sec.update()

    def add(self, cmd):
        prev_valid = self.valid
        prev_design = self.design
        if not prev_valid:
            warnings.warn(
                "This module '{0}' does not currently have valid design string!".format(
                    self.name
                )
            )
        retval = self._set_design(self.design + cmd)
        if not self.valid and prev_valid:
            self.design = prev_design
            raise RuntimeError(
                "Addition to design string '{0}' Makes the design invalid".format(cmd)
            )
        return retval

    def _set_design(self, newdesign):
        self.sec.setDesign(newdesign)
        self.refresh()

    def copyfrom(self, src):
        self.name = src.name
        self.design = src.design
        self.input_switch = src.input_switch
        self.output_switch = src.output_switch
        self.ramp = src.ramp
        self.tolerance = src.tolerance
        self.timeout = src.timeout

    def ZPK_add(
        self,
        filt,
        plane="f",
        annotate=True,
        F_gain_Hz=None,
        scale_gain=True,
    ):
        assert plane in ["s", "f", "n"]

        zpk = TFmath.ANY2ZPKCalc(filt)
        Z, P, K = zpk

        if F_gain_Hz is not None:
            G = 1
            idx = 0
            R = 1
            while True:
                if idx > len(Z) and idx > len(P):
                    break
                if idx < len(Z):
                    R *= Z[idx]
                    G = G * (1j * F_gain_Hz - Z[idx])
                if idx < len(P):
                    R /= P[idx]
                    G = G / (1j * F_gain_Hz - P[idx])
                idx += 1
            if plane == "f":
                K = abs(K / G * (2 * np.pi) ** (len(P) - len(Z)))
            elif plane == "s":
                K = abs(K / G)
            elif plane == "n":
                raise RuntimeError("F_gain_Hz not supported for normalized plane")

        zpk_str = conversion.filter2fotonZPK(
            (Z, P, K),
            plane=plane,
            annotate_pairs=annotate,
            scale_gain=scale_gain,
        )

        return self.add(zpk_str)

    def ZPK_set(self, ZPK, check_response=True, F_response_Hz=None, **kwargs):
        # TODO, add design confirmations that the response is preserved
        self.design = ""
        return self.ZPK_add(ZPK, **kwargs)

    def ZPKz_add(self, ZPK, F_nyquist_Hz=None):
        if F_nyquist_Hz is not None:
            assert F_nyquist_Hz * 2 == self.rate

        def matlab_complex_str(num):
            if np.imag(num) == 0.0:
                return "{0}".format(num)
            else:
                return "{0} + {1}*i".format(np.real(num), np.imag(num))

        Z, P, K = ZPK

        ZPK_template = "zroots([{zeros}],[{poles}],\n{gain})"
        pole_list = []
        zero_list = []
        for pole in P:
            pole_list.append(matlab_complex_str(pole))
        for zero in Z:
            pole_list.append(matlab_complex_str(zero))

        zpk_str = ZPK_template.format(
            zeros="\n\t" + ";\n\t".join(zero_list) + "\n",
            poles="\n\t" + ";\n\t".join(pole_list) + "\n",
            gain=K,
        )
        return self.add(zpk_str)

    def zroots_add(self, ZPK, **kwargs):
        return self.ZPKz_add(ZPK, **kwargs)

    def zroots_set(self, ZPK, **kwargs):
        return self.ZPKz_set(ZPK, **kwargs)

    def gain_add(self, gain, usedb=False):
        """
        @param g gain.
        @param format scalar or dB
        @return true if successful
        """
        if usedb:
            return self.add('gain({0}, "db")\n'.format(gain))
        else:
            return self.add('gain({0}, "scalar")\n'.format(gain))

    def gain_set(self, *args, **kwargs):
        self.design = ""
        return self.gain_add(*args, **kwargs)

    gain_set.__doc__ = gain_add.__doc__

    def butter_add(self, ftype="LowPass", order=None, F1_Hz=None, F2_Hz=None):
        """
        @param type Filter type
        @param order Filter order
        @param f1 Pass band edge (Hz)
        @param f2 Another pass band edge (Hz)
        @return true if successful
        """
        kw = dict(
            ftype=ftype,
            order=order,
            F1_Hz=F1_Hz,
            F2_Hz=F2_Hz,
        )
        if check_ftype(ftype=ftype, F1_Hz=F1_Hz, F2_Hz=F2_Hz):
            return self.add(
                'butter("{ftype}", {order}, {F1_Hz}, {F2_Hz})\n'.format(**kw)
            )
        else:
            return self.add('butter("{ftype}", {order}, {F1_Hz})\n'.format(**kw))

    def butter_set(self, *args, **kwargs):
        self.design = ""
        return self.butter_add(*args, **kwargs)

    butter_set.__doc__ = butter_add.__doc__

    def resgain_add(self, f0=None, Q=None, height=None):
        """
        @param f0 Center frequency.
        @param Q Quality factor ( Q = (Center freq)/(Width) ).
        @param height Height of the peak (dB).
        """
        return self.add(
            "resgain({f0}, {Q}, {height})\n".format(
                f0=f0,
                Q=Q,
                height=height,
            )
        )

    def resgain_set(self, *args, **kwargs):
        self.design = ""
        return self.resgain_add(*args, **kwargs)

    resgain_set.__doc__ = resgain_add.__doc__

    def notch_add(self, f0=None, Q=None, depth=None):
        """
        @param f0 Center frequency.
        @param Q Quality factor ( Q = (Center freq)/(Width) ).
        @param depth depth of the peak (dB).
        """
        return self.add(
            "notch({f0}, {Q}, {height})\n".format(
                f0=f0,
                Q=Q,
                depth=depth,
            )
        )

    def notch_set(self, *args, **kwargs):
        self.design = ""
        return self.notch_add(*args, **kwargs)

    notch_set.__doc__ = notch_add.__doc__

    def elliptic_add(
        self, ftype="LowPass", order=None, rp=None, F1_Hz=None, F2_Hz=None
    ):
        """ """
        check_ftype(ftype=ftype, F1_Hz=F1_Hz, F2_Hz=F2_Hz)
        kw = dict(
            ftype=ftype,
            order=order,
            F1_Hz=F1_Hz,
            F2_Hz=F2_Hz,
        )
        if check_ftype(ftype=ftype, F1_Hz=F1_Hz, F2_Hz=F2_Hz):
            return self.add(
                'elliptic("{ftype}", {order}, {F1_Hz}, {F2_Hz})\n'.format(**kw)
            )
        else:
            return self.add('elliptic("{ftype}", {order}, {F1_Hz})\n'.format(**kw))

    def elliptic_set(self, *args, **kwargs):
        self.design = ""
        return self.elliptic_add(*args, **kwargs)

    elliptic_set.__doc__ = elliptic_add.__doc__

    def cheby1_add(self, ftype="LowPass", order=None, rp=None, F1_Hz=None, F2_Hz=None):
        """
        @param type Filter type
        @param order Filter order
        @param as Stop band attenuation (dB)
        @param f1 Pass band edge (Hz)
        @param f2 Another pass band edge (Hz)
        """
        kw = dict(
            ftype=ftype,
            order=order,
            F1_Hz=F1_Hz,
            F2_Hz=F2_Hz,
        )
        if check_ftype(ftype=ftype, F1_Hz=F1_Hz, F2_Hz=F2_Hz):
            return self.add(
                'cheby1("{ftype}", {order}, {F1_Hz}, {F2_Hz})\n'.format(**kw)
            )
        else:
            return self.add('cheby2("{ftype}", {order}, {F1_Hz})\n'.format(**kw))

    def cheby1_set(self, *args, **kwargs):
        self.design = ""
        return self.cheby1_add(*args, **kwargs)

    cheby1_set.__doc__ = cheby1_add.__doc__

    def cheby2_add(self, ftype="LowPass", order=None, rp=None, F1_Hz=None, F2_Hz=None):
        """
        @param type Filter type
        @param order Filter order
        @param f1 Pass band edge (Hz)
        @param f2 Another pass band edge (Hz)
        @return true if successful
        """
        kw = dict(
            ftype=ftype,
            order=order,
            F1_Hz=F1_Hz,
            F2_Hz=F2_Hz,
        )
        if check_ftype(ftype=ftype, F1_Hz=F1_Hz, F2_Hz=F2_Hz):
            return self.add(
                'cheby2("{ftype}", {order}, {F1_Hz}, {F2_Hz})\n'.format(**kw)
            )
        else:
            return self.add('cheby2("{ftype}", {order}, {F1_Hz})\n'.format(**kw))

    def cheby2_set(self, *args, **kwargs):
        self.design = ""
        return self.cheby2_add(*args, **kwargs)

    cheby2_set.__doc__ = cheby2_add.__doc__


def iir2zpkstr(filt, plane="s", prewarp=True):
    """Returns the zeros and poles of an IIR filter.
    The returned string has the format "zpk(...)", if the plane is
    "s", "n" or "f". It is of the form "rpoly(...)", if the plane is
    "p".
    """
    try:
        arg = filt.filt.get()
    except AttributeError:
        arg = filt.get()
    zpk = ROOT.string()
    ROOT.iir2zpk(arg, zpk, plane, prewarp)
    return zpk


def iir2zpk(filt, plane="s", prewarp=True):
    """Returns the zeros and poles of an IIR filter.
    The returned string has the format "zpk(...)", if the plane is
    "s", "n" or "f". It is of the form "rpoly(...)", if the plane is
    "p".
    """
    try:
        arg = filt.filter_raw.get()
    except AttributeError:
        arg = filt.get()
    assert plane in ["s", "f", "n"]

    nsos = ROOT.iirsoscount(arg)
    Nz = ctypes.c_int(0)
    Np = ctypes.c_int(0)
    arr_z = np.empty(4 * nsos, np.double)
    arr_p = np.empty(4 * nsos, np.double)
    gain = ROOT.Double(0)
    ROOT.iir2zpk(
        arg,
        Nz,
        cppyy.bind_object(arr_z, TcplxD),
        Np,
        cppyy.bind_object(arr_p, TcplxD),
        gain,
        plane,
        prewarp,
    )
    Np = Np.value
    Nz = Nz.value
    poles = arr_p[: 2 * Np : 2] + 1j * arr_p[1 : 2 * Np : 2]
    zeros = arr_z[: 2 * Nz : 2] + 1j * arr_z[1 : 2 * Nz : 2]
    return zeros, poles, gain


def iir2xfer(filt, F_Hz):
    """Returns the zeros and poles of an IIR filter.
    The returned string has the format "zpk(...)", if the plane is
    "s", "n" or "f". It is of the form "rpoly(...)", if the plane is
    "p".
    """
    try:
        arg = filt.filter_raw.get()
    except AttributeError:
        arg = filt.get()
    F_Hz = np.asarray(F_Hz, dtype=np.float32, order="C")
    arr_H = np.empty(len(F_Hz) * 2, dtype=np.float32, order="C")

    arg.Xfer(
        cppyy.bind_object(arr_H, TcplxF),
        F_Hz,
        len(F_Hz),
    )

    return arr_H[::2] + 1j * arr_H[1::2]


def iir2sos(filt, format="s"):
    """Returns the a's and b's of an IIR filter.
    The returned a's and b's are grouped in zecond order sections.
    The first returned coeffcient is the overall gain. The following
    coeffcients come in groups of 4. If the format is 's' (standard),
    the order is b1, b2, a1, b2. If the format is 'o' (online), the
    order is a1, a2, b1, b2. The returned length is always of the
    for nba = 1 + 4 * (number of second order sections). A second order
    section is defined as:
    $$H(z)=\\frac{b0+b1 z^{-1}+b2 z^{-2}}{1+a1 z^{-1}+a2 z^{-2}}$$
    """
    try:
        arg = filt.filter_raw.get()
    except AttributeError:
        arg = filt.get()
    nsos = ROOT.iirsoscount(arg)
    nba = ctypes.c_int(0)
    ba = array("d", range(1 + 4 * nsos))
    ROOT.iir2z(arg, nba, ba, format)
    return [ba[n] for n in range(nba)]


def iir2zroots(filt):
    """Returns the a's and b's of an IIR filter.
    The returned a's and b's are grouped in zecond order sections.
    The first returned coeffcient is the overall gain. The following
    coeffcients come in groups of 4. If the format is 's' (standard),
    the order is b1, b2, a1, b2. If the format is 'o' (online), the
    order is a1, a2, b1, b2. The returned length is always of the
    for nba = 1 + 4 * (number of second order sections). A second order
    section is defined as:
    $$H(z)=\\frac{b0+b1 z^{-1}+b2 z^{-2}}{1+a1 z^{-1}+a2 z^{-2}}$$
    """
    try:
        arg = filt.filter_raw.get()
    except AttributeError:
        arg = filt.get()
    nsos = ROOT.iirsoscount(arg)
    nba = ctypes.c_int(0)
    ba = array("d", range(1 + 4 * nsos))
    ROOT.iir2z(arg, nba, ba, "s")
    Nz = ctypes.c_int(0)
    Np = ctypes.c_int(0)
    arr_z = np.empty(4 * nsos, np.double)
    arr_p = np.empty(4 * nsos, np.double)
    gain = ROOT.Double(0)
    ROOT.z2z(
        nba,
        ba,
        Nz,
        cppyy.bind_object(arr_z, TcplxD),
        Np,
        cppyy.bind_object(arr_p, TcplxD),
        gain,
        "s",
    )
    Np = Np.value
    Nz = Nz.value
    poles = arr_p[: 2 * Np : 2] + 1j * arr_p[1 : 2 * Np : 2]
    zeros = arr_z[: 2 * Nz : 2] + 1j * arr_z[1 : 2 * Nz : 2]
    return zeros, poles, gain
