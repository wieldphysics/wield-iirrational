from IIRrational.utilities.ipynb_lazy import *
from IIRrational.utilities.ipynb.sympy import *
from IIRrational.TFmath import order_reduce
import numpy as np
import IIRrational
import scipy
from scipy import signal
from os import path


#z = eig([A B;C D],diag([ones(1,n) 0]);



def test_statespace_fit(tpath):
    cpath = path.split(__file__)[0]
    ss = IIRrational.load(path.join(cpath, 'HSTS.mat'))['HSTS']
    F_Hz = np.logspace(-1, +1, 2000)
    A, B, C, D = ss['A'], ss['B'], ss['C'], ss['D']

    idx_in = 0
    idx_out = 2

    fit = IIRrational.v2.ss2filter(
        A, B, C, D,
        F_Hz = F_Hz,
        idx_in = idx_in,
        idx_out = idx_out
    )
    fit.choose(10)

    axB = fit.investigate_fit_plot()
    axB.save(path.join(tpath, 'plot.pdf'))


