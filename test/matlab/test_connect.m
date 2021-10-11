% Author : Lee McCuller
clear classes
addpath(char(py.IIRrational.matlabpath()));

%IIRrational = msurrogate.PySurrogate();
%IIRrational.connect_subprocess('IIRrational.matlab')

iir = IIRrational.surrogate();
dat = iir.testing.IIRrational_data({'simple1'})
