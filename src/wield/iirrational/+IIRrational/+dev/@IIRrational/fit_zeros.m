function retS = fit_zeros(varargin)
  p = inputParser ();
  p.FunctionName = "fit_zeros";
  p.addRequired('F_Hz', @(v)~ischar(v))
  p.addRequired('data', @(v)~ischar(v))
  p.addParamValue('W', 1)
  p.addParamValue('ZPK', {[],[],1})
  p.addParamValue('F_nyquist', 1)
  p.addParamValue('nzeros', 30)
  p.addParamValue('npoles', 30)
  p.addParamValue('compute_residuals', 1)
  p.addParamValue('stabilize', 1)
  p.addParamValue('mindelay',  0)
  p.parse(varargin{:});
  args = p.Results;
  retS = args;
    #X, poly, ZPK = XpolyZPK(F_Hz, F_nyquist, ZPK)
  zeros_init = args.ZPK{1};
  poles_init = args.ZPK{2};
  gain_init = args.ZPK{3};

  X = exp(1j*pi*args.F_Hz / args.F_nyquist);
  V_b = vander(X, args.nzeros + 1);
  a_init = poly(poles_init);
  ah = polyval(a_init, X);
  retS = ah;

    #A_b = np.einsum('ij,i->ij', X_b,  W / (data * ah))
    #b_fit, res, rank, s = scipy.linalg.lstsq(
    #    np.vstack([A_b.real, A_b.imag]),
    #    np.hstack([W.real, W.imag]),
    #)
    ##not needed when projecting to the purely real problem
    ##poly.coeffclean(b_fit)
    #gain = b_fit[0]
    #zeros = poly.roots(b_fit)

    #ZPK = poly.root_transform(
    #    (zeros, poles_init, gain),
    #    stabilize = False,
    #    mindelay  = mindelay,
    #)

    #retB = QuickBunch(
    #    ZPK         = ZPK,
    #    data        = data,
    #    F_Hz       = F_Hz,
    #    W       = W,
    #    F_nyquist_Hz   = F_nyquist,
    #)

    #if compute_residuals:
    #    bh = poly.val(X, b_fit)
    #    abh = bh / ah
    #    retB.data_fit    = abh
    #    debias_reweight = 1/(.001 + W**2)
    #    retB.resP = W * (abh/data - 1)
    #    retB.resZ = W * (data/abh - 1) * debias_reweight
    #    retB.res = np.sum((abs(retB.resP)**2 + abs(retB.resZ)**2) / (1 + debias_reweight)) / (2 * len(data))
    #return retB
endfunction



#def fit_poles_mod_zeros(
#        F_Hz,
#        data,
#        W         = 1,
#        Wmod      = None,
#        npoles    = None,
#        nzeros    = None,
#        ZPK       = ((), (), 1),
#        F_nyquist_Hz = None,
#        max_size  = None,
#        stabilize = True,
#        **kwargs
#):
#    F_Hz, data, W = np.broadcast_arrays(F_Hz, data, W)
#
#    if Wmod is None:
#        Wmod = W
#    else:
#        W, Wmod = np.broadcast_arrays(W, Wmod)
#
#    X, poly, ZPK = XpolyZPK(F_Hz, F_nyquist, ZPK)
#    zeros_init, poles_init, gain_init = ZPK
#
#    X_b = poly.vander(X, nzeros + 1)
#    X_a = poly.vander(X, npoles + 1)
#    bh = poly.valfromroots(X, zeros_init)
#    ah = poly.valfromroots(X, poles_init)
#
#    A_b = np.einsum('ij,i->ij', X_b,  Wmod / (data * ah))
#    #remove the effect of linear (delay) phasing
#    #A_b = np.hstack([(1j*F_Hz).reshape(-1, 1), A_b])
#
#    q, r = np.linalg.qr(A_b)
#
#    A_a = np.einsum('ij,i->ij', X_a,  W * data / bh)
#    A_a = A_a - np.einsum('ij,jk->ik', q, np.einsum('ij,ik->jk', q.conjugate(), A_a))
#
#    S, V = SVD_SV(
#        np.vstack([A_a.real, A_a.imag]),
#        n_smallest  = 3,
#        overwrite_a = True,
#    )
#    #print("POLES SVD: ", S[-4:])
#    a_fit = V.T[:, -1]
#    #not need with fully real problem
#    #poly.coeffclean(a_fit)
#
#    gain = 1/a_fit[0]
#
#    poles = poly.roots(a_fit)
#    ZPK = poly.root_transform(
#        (zeros_init, poles, gain),
#        stabilize = stabilize
#    )
#
#    retB = QuickBunch(
#        ZPK         = ZPK,
#        data        = data,
#        F_Hz       = F_Hz,
#        W       = W,
#        F_nyquist_Hz   = F_nyquist,
#    )
#
#    return retB
#
#
#def fit_zeros_mod_poles(
#        F_Hz,
#        data,
#        W         = 1,
#        Wmod      = None,
#        npoles    = None,
#        nzeros    = None,
#        ZPK       = ((), (), 1),
#        F_nyquist_Hz = None,
#        max_size  = None,
#        mindelay  = True,
#        **kwargs
#):
#    F_Hz, data, W = np.broadcast_arrays(F_Hz, data, W)
#
#    if Wmod is None:
#        Wmod = W
#    else:
#        W, Wmod = np.broadcast_arrays(W, Wmod)
#
#    X, poly, ZPK = XpolyZPK(F_Hz, F_nyquist, ZPK)
#    zeros_init, poles_init, gain_init = ZPK
#
#    X_b = poly.vander(X, nzeros + 1)
#    X_a = poly.vander(X, npoles + 1)
#    bh = poly.valfromroots(X, zeros_init)
#    ah = poly.valfromroots(X, poles_init)
#
#    A_a = np.einsum('ij,i->ij', X_a,  Wmod / (data / bh))
#    #remove the effect of linear (delay) phasing
#    #A_a = np.hstack([(1j*F_Hz).reshape(-1, 1), A_a])
#
#    q, r = np.linalg.qr(A_a)
#
#    A_b = np.einsum('ij,i->ij', X_b,  W / data / ah)
#    A_b = A_b - np.einsum('ij,jk->ik', q, np.einsum('ij,ik->jk', q.conjugate(), A_b))
#    S, V = SVD_SV(
#        np.vstack([A_a.real, A_a.imag]),
#        n_smallest  = 3,
#        overwrite_a = True,
#    )
#    #print("ZEROS SVD: ", S[-4:])
#    b_fit = V.T[:, -1]
#    #not needed with fully real problem
#    #poly.coeffclean(b_fit)
#
#    gain = 1/b_fit[0]
#
#    zeros = poly.roots(b_fit)
#    ZPK = poly.root_transform(
#        (zeros, poles_init, gain),
#        stabilize = False,
#        mindelay  = mindelay,
#    )
#
#    retB = QuickBunch(
#        ZPK       = ZPK,
#        data      = data,
#        F_Hz      = F_Hz,
#        W         = W,
#        F_nyquist_Hz = F_nyquist,
#    )
#
#    return retB
#
#def fit_poles(
#    F_Hz,
#    data,
#    npoles,
#    W                 = 1,
#    nzeros            = None,
#    F_nyquist_Hz         = None,
#    ZPK               = ((), (), 1),
#    compute_residuals = True,
#    stabilize         = False,
#    **kwargs
#):
#    F_Hz, data, W = np.broadcast_arrays(F_Hz, data, W)
#
#    X, poly, ZPK = XpolyZPK(F_Hz, F_nyquist, ZPK)
#    zeros_init, poles_init, gain_init = ZPK
#
#    X_a = poly.vander(X, npoles + 1)
#    bh = poly.valfromroots(X, zeros_init)
#
#    A_a = np.einsum('ij,i->ij', X_a,  W * data / bh)
#    #solve the problem with purely real taps
#    a_fit, res, rank, s = scipy.linalg.lstsq(
#        np.vstack([A_a.real, A_a.imag]),
#        np.hstack([W.real, W.imag]),
#    )
#    #not needed when projecting to the purely real problem
#    #poly.coeffclean(a_fit)
#    gain = 1/a_fit[0]
#    poles = poly.roots(a_fit)
#
#    ZPK = poly.root_transform(
#        (zeros_init, poles, gain),
#        stabilize = stabilize,
#    )
#
#    retB = QuickBunch(
#        ZPK       = ZPK,
#        data      = data,
#        F_Hz      = F_Hz,
#        W         = W,
#        F_nyquist_Hz = F_nyquist,
#    )
#
#    if compute_residuals:
#        ah = poly.val(X, a_fit)
#        abh = bh / ah
#        retB.data_fit    = abh
#        debias_reweight = 1/(.001 + W**2)
#        retB.resP = W * (abh/data - 1)
#        retB.resZ = W * (data/abh - 1) * debias_reweight
#        retB.res = np.sum((abs(retB.resP)**2 + abs(retB.resZ)**2) / (1 + debias_reweight)) / (2 * len(data))
#    return retB

