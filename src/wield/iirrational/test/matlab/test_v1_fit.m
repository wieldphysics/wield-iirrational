% Author : Lee McCuller
clear classes
addpath(char(py.wield.iirrational.matlabpath()));

iir = wield.iirrational.surrogate();
dat = iir.testing.IIRrational_data({'simple1'});

%now convert to native
kw = struct();
kw.data = dat.data;
kw.F_Hz = dat.F_Hz;
kw.F_nyquist_Hz = dat.F_nyquist_Hz;

fit_out = iir.v1.data2filter(...
  kw, ...
  'SNR', dat.SNR    ...
                            );

%axB = iir.plots.plot_fitter_flag({fit_out.fitter});
axB = iir.plots.plot_fitter_flag({fit_out.fitter}, 'fname', 'out-simple1/test.png');
