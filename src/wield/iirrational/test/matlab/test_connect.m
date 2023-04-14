% Author : Lee McCuller
clear classes
addpath(char(py.wavestate.iirrational.matlabpath()));

%wavestate.iirrational = msurrogate.PySurrogate();
%wavestate.iirrational.connect_subprocess('wavestate.iirrational.matlab')

iir = wavestate.iirrational.surrogate();
dat = iir.testing.IIRrational_data({'simple1'})
