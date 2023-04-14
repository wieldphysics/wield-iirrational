% Log plot of magnitude and phase of a complex vector.
%
% zplotlog(f, h)

function handle = zplotmag(f, h, varargin)

  % phase
  subplot(2, 1, 2);
  semilogx(f, 180 * angle(h) / pi, varargin{:});
  ylabel('phase [deg]');
  grid();

  % magnitude (done second so that it is "selected")
  subplot(2, 1, 1);
  handle = semilogx(f, abs(h), varargin{:});
  ylabel('mag');
  grid();


