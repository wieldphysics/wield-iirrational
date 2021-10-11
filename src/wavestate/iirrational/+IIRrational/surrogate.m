%
%
function iirsur = surrogate(optional_cookie)
  persistent sur;

  if ~exist('msurrogate')
      addpath(char(py.msurrogate.matlabpath()));
  end

  if isempty(sur)
    sur = msurrogate.PySurrogate();
  end

  iirsur = sur;

  switch nargin
  case 1
    if ~isempty(optional_cookie)
      iirsur.connect_cookie(optional_cookie);
    end
  case 0
    kw = struct();
    kw.env.LD_LIBRARY_PATH = '';
    if ~iirsur.attached()
      iirsur.connect_subprocess('IIRrational.matlab', kw);
    end
  end
end
