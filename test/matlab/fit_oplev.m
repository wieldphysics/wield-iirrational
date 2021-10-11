%% 


%% Start up iirrational
addpath(char(py.iirrational.matlabpath()))
iir = iirrational.surrogate()

%%

ZPKs = struct()

%% now load data

AP = open('ArnaudOplevPlant.mat');

name = 'ix'
fname = 'ix2.pdf'
subdata = AP.ap.(name);

ff = subdata.ff;
SNR = 100 * (ff < 3) + 3 * (ff < 50);

select = 1:2:900;
kw = struct();
kw.data = subdata.plant(select).';
kw.F_Hz = subdata.ff(select).';
kw.SNR = SNR(select).';
kw.F_nyquist_Hz = 16384/2.
kw.alt_res = true;
kw.order_initial = 15;

%fit it out
out = iir.v1.data2filter(kw)

select = 1:5:1000;
kw = struct();
kw.data = subdata.plant(select).';
kw.F_Hz = subdata.ff(select).';
kw.SNR = SNR(select).';
kw.F_nyquist_Hz = 16384/2.
kw.ZPK = out.fitter.ZPK;
%kw.alt_res = true;
%out = iir.v1.data2filter(kw)


iir.plots.plot_fitter_flag({out.fitter}, 'fname', fname)


disp('Z')
out.fitter.ZPK{1}
disp('P')
out.fitter.ZPK{2}
disp('K')
out.fitter.ZPK{3}
ZPK.(name) = out.fitter.ZPK


