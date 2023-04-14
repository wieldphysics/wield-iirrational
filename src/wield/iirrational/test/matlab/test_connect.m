% Author : Lee McCuller
clear classes
addpath(char(py.wield.iirrational.matlabpath()));

%wield.iirrational = msurrogate.PySurrogate();
%wield.iirrational.connect_subprocess('wield.iirrational.matlab')

iir = wield.iirrational.surrogate();
dat = iir.testing.IIRrational_data({'simple1'})
