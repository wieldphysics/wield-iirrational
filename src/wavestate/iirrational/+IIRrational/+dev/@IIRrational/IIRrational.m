%% @functiontypefn  {Function File} {} polynomial ()
%% @functiontypefnx {Function File} {} polynomial (@var{a})
%% Create a polynomial object representing the polynomial
%%
%% @example
%% @end example
%%
%% @noindent
%% @end functiontypefn
classdef RationalDiscFilter < DepStruct
  %properties
  %end
  methods
    function self = RationalDiscFilter(varargin)
      self = self@DepStruct(varargin{:});
    end

  end
  methods (Access = public)
    function build(self, varargin)
      p = inputParser();
      p.FunctionName = 'build';
      p.addOptional('F_Hz', 1, @(v)~ischar(v));
      p.addOptional('data', 1, @(v)~ischar(v));
      p.addParamValue('W', 1);
      p.addParamValue('ZPK', {[],[],1});
      p.addParamValue('F_nyquist', 1);
      p.addParamValue('nzeros', 30);
      p.addParamValue('npoles', 30);
      p.parse(varargin{:});
      args = p.Results;

      self.add_generator('a_vec', @(self) self.gain * poly(self.poles))

      function val = a_vec(self, val)
        self.dependencies('poles', 'gain');
        if self.get_raw('stabilize');
            %TODO need to fix
            %self.gain  =
            self.poles = self.poly.roots(val);
            val = self.TAG_NO_SET;
          end
      end
      self.add_setter(@a_vec);

      self.add_generator('b_vec', @(self) self.gain * poly(self.poles))

      function val = b_vec(self, val)
        self.dependencies('zeros', 'gain');
        if self.get_raw('mindelay');
              %actually go straight to poles since that is where they are stabilized
              %TODO need to fix
              %self.gain  =
          self.zeros = self.poly.roots(val);
          val = self.TAG_NO_SET;
        end
      end
      self.add_setter(@b_vec);

      self.add_generator('zeros', @(self) poly(self.b_vec));

      function val = zeros(self, val)
        if self.get_raw('mindelay');
          val = self.root_transform(val);
        end
        self.dependencies('b_vec');
      end
      self.add_setter(@zeros);

      self.add_generator('poles', @(self) poly(self.a_vec));

      function val = poles(self, val)
        if self.get_raw('stabilize');
          val = self.root_transform(val);
        end
        self.dependencies('a_vec');
      end
      self.add_setter(@poles);

      self.add_generator('gain', @(self) self.b_vec(1) ./ self.a_vec(1));

      self.add_generator('V_b', @(self) vander(self.X_grid, self.nzeros + 1));
      self.add_generator('V_a', @(self) vander(self.X_grid, self.npoles + 1));

      self.add_generator('h_a', @(self) polyval(self.X_grid, self.a_vec));
      self.add_generator('h_b', @(self) polyval(self.X_grid, self.b_vec));

      self.add_generator('data_fit', @(self) self.h_b / self.h_a);

      function retB = residuals(self)
        debias_reweight = 1/(.001 + self.W^2);
        retB.resP = self.W .* (self.data_fit/self.data - 1);
        retB.resZ = self.W .* (self.data/self.data_fit - 1);
        retB.total = sum(...
          (abs(retB.resP)^2 + abs(retB.resZ * debias_reweight)^2) / (1 + debias_reweight)...
                        ) / (2 * len(self.data));
      end
      self.add_generator(@residuals);

      self.add_generator('A_b', @(self) self.V_b * (self.W ./ (self.data .* self.h_a)));
      self.add_generator('A_a', @(self) self.V_a * (self.W .* (self.data ./ self.h_b)));

      self.add_generator('X_grid', @(self) exp(1j * pi * self.F_Hz / self.F_nyquist));

      %needs more error checking!
      self.my.mindelay  = false;
      self.my.stabilize = false;
      self.my.F_Hz      = args.F_Hz;
      self.my.data      = args.data;
      self.my.W         = args.W;
      self.my.F_nyquist_Hz = args.F_nyquist;
      self.my.zeros     = args.ZPK{1};
      self.my.poles     = args.ZPK{2};
      self.my.gain      = args.ZPK{3};
    end

    function roots = root_transform(self, roots)
      select = (abs(roots) > 1);
      roots(select) = 1 / conjugate(roots(select))
    end

    function fit_poles_mod_zeros(self, varargin)
      %remove the effect of linear (delay) phasing
      A_b = self.A_b;
      %A_b = np.hstack([(1j*F_Hz).reshape(-1, 1), A_b])
      q, r = qr(A_b);

      A_a = self.A_a;
      A_a = A_a - np.einsum('ij,jk->ik', q, np.einsum('ij,ik->jk', q.conjugate(), A_a));

      [U, S, V] = SVD([A_a.real, A_a.imag], 'econ');
      %print("POLES SVD: ", S[-4:])
      self.a_vec = V.T(:, end);
    end

    function fit_zeros_mod_poles(self, varargin)
        A_a = self.A_a
        %remove the effect of linear (delay) phasing
        %A_a = np.hstack([(1j*F_Hz).reshape(-1, 1), A_a])

        q, r = qr(A_a)

        A_b = self.A_b
        A_b = A_b - np.einsum('ij,jk->ik', q, np.einsum('ij,ik->jk', q.conjugate(), A_b))
        S, V = SVD([A_a.real, A_a.imag], 'econ');
        %print("ZEROS SVD: ", S[-4:])
        self.b_vec = V.T(:, end)
    end

    function fit_poles(self, varargin)
      %print(self.a_vec, self.b_vec)
      A_a = self.A_a;
      %solve the problem with purely real taps
      opts.RECT=true;
      a_fit, res, rank, s = linsolve([A_a.real, A_a.imag], [self.W.real, self.W.imag], opts);
      self.a_vec = a_fit;
    end

    function fit_zeros(self, varargin)
      A_b = self.A_b;
      opts.RECT=true;
      b_fit, res, rank, s = linsolve([A_b.real, A_b.imag], [self.W.real, self.W.imag], opts);
      self.b_vec = b_fit;
    end

    function fit_pzpz(self, varargin)
      %max_size          = None,
      %collect_all       = False,
      %zeros_first       = False,
      %stabilize         = False,
      %stabilize_svd     = True,
      %mindelay          = False,
      %mindelay_svd      = False,
      %n_svd             = 1,
      %n_iter            = 10,
      pre.mindelay  = self.mindelay
      pre.stabilize = self.stabilize

      collection = []
      if not(zeros_first)
          if n_svd >= 1
              fitA = self.fit_poles_mod_zeros
          else
              fitA = self.fit_poles
          end
          if n_svd >= 2
              fitB = self.fit_zeros_mod_poles
          else
              fitB = self.fit_zeros
          end
          fitC = self.fit_poles
          fitD = self.fit_zeros
      else
          if n_svd >= 1
              fitA = self.fit_zeros_mod_poles
          else
              fitA = self.fit_zeros
          end
          if n_svd >= 2
              fitB = self.fit_poles_mod_zeros
          else
              fitB = self.fit_poles
          end
          fitC = self.fit_zeros
          fitD = self.fit_poles
      end

      self.mindelay  = mindelay_svd
      self.stabilize = stabilize_svd
      fitA('max_size', max_size)
      fitB('max_size', max_size)
      self.mindelay  = mindelay
      self.stabilize = stabilize

      for i = 1:n_iter
          fitC()
          fitD()
          %if collect_all
              %collection.append(self.copy())
      end

      if n_iter > 0
          fitC()
          fitD()
      end

      self.stabilize = pre.stabilize
      self.mindelay  = pre.mindelay
      %if collect_all
          %collection.append(self.copy())
    end

    function fit_pz(self, varargin)
      %n_iter = 0,
      %n_svd = 1,
      self.fit_pzpz(varargin{:});
    end
  end
end
