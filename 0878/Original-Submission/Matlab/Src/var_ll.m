% VAR_LL  Vector-autoregressive likelihood and optionally its derivative
%
%   [LL, OK] = VAR_LL(X, A, Sig) sets LL to the exact log-likelihood function
%   value for n observations of the r-variate zero-mean VAR time series:
%
%            x(t) = A1¸x(t-1) + ... + Ap¸x(t-p) + eps(t)
%
%   where x(t) is r-dimensional and eps(t) is r-variate normal with mean 0 and
%   covariance Sig. Sig and the Ai's are rÚr matrices. A should contain the r Ú
%   r¸p matrix [A1 A2...Ap] and X should have x(t) in its t-th column for t =
%   1,...,n. To use this likelihood routine to fit a non-zero-mean time-series,
%   the mean-vector of the observations may be subtracted from each x(t). OK =
%   true indicates success, but is false if the model is non-stationary.
%
%   [LL, OK] = VAR_LL(X, A, Sig, mu, miss) is used when there are missing
%   values. In that case it is appropriate to let the mean-vector of the series,
%   mu, be a parameter, and use the model:
%
%            x(t) - mu = A1¸(x(t-1) - mu) + ... + Ap¸(x(t-p) - mu) + eps(t)
%
%   The rÚn logical array miss should be true in missing locations. If missing
%   values are NaN in X one can call [LL, OK] = VAR_LL(X, A, Sig, mu, isnan(X)).
%   An empty mu is equivalent to a zero mu.
%
%   [LL, OK, EPS] = VAR_LL(..., 'res') returns the maximum likelihood estimate
%   of the residuals in EPS and [LL, OK, EPS, XM] = VAR_LL(...,'res_miss') finds
%   in addition the maximum likelihood estimate of the missing values in XM.
%
%   [LL, OK, LLD] = VAR_LL(X, A, Sig) and [LL, OK, LLD] = VAR_LL(X, A, Sig, mu,
%   miss) return the log-likelihood function value in LL and its gradient in
%   LLD. If mu is empty (or not present), LL is as if mu were zero, but LLD is
%   without mu-derivatives.
%
%   [LL,OK,LLD] = VAR_LL(...,J) may be used to speed up the gradient calculation
%   when A and/or Sig depend on a smaller set of independent variables and only
%   the gradient w.r.t. this smaller set is sought. J gives the Jacobian of the
%   variable change, and should have either r^2¸p rows or r^2¸p + r¸(r+1)/2
%   rows. In the latter case both A and Sig depend on a smaller set, but in the
%   former case it is only A that does. Thus the possible variable changes are
%   theta --> A and theta --> [A vech(Sig)]. See also help text of varma_llc.
%
%   The so-called "Cholesky decomposition" method is used, see [1] and [2].
%
%   [1] K Jonasson and SE Ferrando 2006. Efficient likelihood evaluation for
%       VARMA processes with missing values. Report VHI-01-2006, Engineering
%       Research Institute, University of Iceland.
%
%   [2] K Jonasson 2006. Matlab programs for complete and incomplete data exact
%       VARMA likelihood and its gradient. Report VHI-02-2006, Engineering
%       Research Institute, University of Iceland.
%
%   Kristjßn Jnasson, Dept. of Computer Science, University of Iceland, 2006.
%   jonasson@hi.is.

function [ll, ok, varargout] = var_ll(X, A, Sig, varargin)
  if nargin>3 && ischar(varargin{end}), code = varargin{end}; else code = '';end
  MISS = nargin >=5 && ~iscell(varargin{1}) && islogical(varargin{2});
  FINDRES = any(strmatch(code, {'res','res_miss'}, 'exact'));
  FNDGRAD = nargout == 3 && ~FINDRES; % find gradient
  CHNGVAR = FNDGRAD && (MISS && nargin>=6 || ~MISS && nargin>=4);
  LowTri.LT = true;
  [p, q, r, n] = get_dimensions(A, [], Sig, X); % ascertain(p>0);
  nPar = r^2*p + r*(r+1)/2;
  if ~MISS
    if nargin > 3 && islogical(varargin{1}), error('wrong arguments'), end
    MEAN = false;
    mu = zeros(r,1);
    miss = false(r,n);
    if CHNGVAR, J = varargin{1}; end
  else
    [mu, miss] = deal(varargin{1:2});
    mu=mu(:);
    MEAN = ~isempty(mu);
    if ~MEAN, mu = zeros(r,1); end
    if CHNGVAR, J = varargin{3}; end
  end
  [ko,ro,km] = find_missing_info(miss);
  mubar = repmat(mu, n, 1);
  X = X(:);
  xo = X(~miss) - mubar(~miss);
  nObs = length(xo);
  if nObs == n*r, MISS = false; end
  ll = 0;
  varargout = {0};
  PLU = vyw_factorize(A);
  vyw_ok = isempty(PLU) || isempty(PLU{1}) || PLU{1}(1)~=0;
  ok = vyw_ok && is_stationary(A, Sig, PLU);
  if ~ok, if nargout<=1, error('Non-stationary model'); else return, end, end
  %
  if ~FNDGRAD || FINDRES % ONLY FUNCTION VALUE
    S              = vyw_solve           (A, PLU, {Sig});
    [Lu,LSig,jpat] = omega_ar_factor     (S, Sig, miss);
    wo             = lambda_multiply     (A, xo, miss);
    wohat          = omega_ar_trisolve   (Lu, LSig, wo, jpat, ko, 'NoT');
    logdetO        = omega_ar_logdet     (Lu, LSig, jpat);
    if MISS
      Gneg         = find_Gneg           (A, [], {Sig}, n);
      Scol         = S_extend            (A, {Sig}, S, n);
      Acol         = find_Acol           (A, r);
      Laom         = find_lambda_om      (Acol, miss); % Lambda_om
      Laomh        = profile_ar_trisolve (Lu, LSig, Laom, jpat, ko, km,p,p,'N');
      clear Laom
      V            = find_V              ({Sig}, Gneg, Scol, miss);
      Vhat         = profile_ar_trisolve (Lu, LSig, V, jpat, ko, km, p, q, 'N');
      clear V
      Rlam         = profile_ata         (Laomh, ko, km, p, p);
      P            = profile_mult        (Laomh, Vhat, ko, km, p, p, q);
      Lw           = profile_mult        (Laomh, wohat, ko, km, p, p, n);
      if ~FINDRES,
        clear Laomh
      end
      RV           = profile_ata         (Vhat, ko, km, p, q);
      Vw           = profile_mult        (Vhat, wohat, ko, km, p, q, n);
      if ~FINDRES
        clear Vhat
      end
      Sm           = find_Sm             (Scol, miss);
      SmP          = Sm*P;                                            clear P
      R            = -SmP-SmP' + atba_c  (Sm, Rlam, Sm + RV);         clear Rlam
      LR           = cholf               (R);                         clear R
      K            = linsolve            (LR, RV - SmP, LowTri);      clear SmP
      Q            = Sm - RV + K'*K;                                  clear RV
      LQ           = cholf               (Q);                         clear Q
      u            = linsolve            (LR, Vw - Sm*Lw, LowTri);
      v            = linsolve            (LQ, Vw - K'*u, LowTri);
      xSxo         = wohat'*wohat - u'*u + v'*v;
      logdetLR     = sum                 (log(diag(LR)));
      logdetLQ     = sum                 (log(diag(LQ)));
      logdetSm     = 2*sum               (log(diag(cholf(Sm))));
      logdetSo     = logdetO + 2*logdetLR + 2*logdetLQ - 2*logdetSm;
    else
      xSxo         = wohat'*wohat;
      logdetSo     = logdetO;
      LR = []; LQ = []; K = []; u = []; v = []; Laomh = []; Vhat = []; Sm = [];
    end

    ll = -0.5*(nObs*log(2*pi) + logdetSo + xSxo);
    
    if FINDRES && MISS
      [eps,xm]=res_miss_ar(A,Sig,Lu,LSig,wohat,LR,LQ,K,u,v,Laomh,Vhat,Sm,miss);
      varargout = {eps, xm + mubar(miss)};
    elseif FINDRES
      [eps,xm] = res_miss_ar(A,Sig,Lu,LSig,wohat);
      varargout = {eps, xm};
    end

  else  % FIND ALSO GRADIENT
    if MEAN
      xod = xmu_deriv(nPar + r, miss);
    else
      xod = zeros(nObs, 1, nPar);
    end
    [wo, wod]      = lambda_multiply   (A, xo, miss, xod);
    SigS           = mds_set_parmat    (p+1, Sig, p+1);
    Sigd           = der2array         (SigS);
    [C,G,Cd,Gd]    = deal              ({Sig}, {Sig}, Sigd, Sigd);
    S              = vyw_solve         (A, PLU, {Sig});
    RHS            = vyw_deriv_rhs     (A, SigS, S);
    Sd             = vyw_solve         (A, PLU, RHS);
    if MISS
      [Gneg, Gnegd]  = find_Gneg         (A, [], C, n, Cd);
      [Scol, Scold]  = S_extend          (A, G, S, n, Gd, Sd);
      [Acol, Acold]  = find_Acol         (A, r, nPar);
    end
    if CHNGVAR
      wod    = chng_var(wod, J);
      Sigd   = chng_var(Sigd, J);
      Gd     = chng_var(Gd, J);
      Sd     = chng_var(Sd, J);
      if MISS
        Gnegd  = chng_var(Gnegd, J);
        Scold  = chng_var(Scold, J);
        Acold  = chng_var(Acold, J);
      end
    end
    [Lu,LSig,jpat,Lud,LSigd] = omega_ar_factor(S, Sig, miss, Sd, Sigd{1});
    if MISS
      [Laom, Laomd]    = find_lambda_om      (Acol, miss, Acold); % Lambda_om
      [Laomh,Laomhd]   = profile_ar_trisolve (Lu, LSig, Laom, jpat, ko, km, p...
        ,                                            p, 'N', Lud, LSigd, Laomd);
      clear Laom Laomd
      [V, Vd]          = find_V              (G, Gneg,Scol,miss,Gd,Gnegd,Scold);
      [Vhat, Vhatd]    = profile_ar_trisolve (Lu, LSig, V, jpat, ko, km, p, q...
        ,                                                  'N', Lud, LSigd, Vd);
      clear V Vd
    end
    if MEAN, [Lud,LSigd] = add_mu_deriv  (r, Lud, LSigd); end
    [woh, wohd]     = omega_ar_trisolve  (Lu,LSig,wo,jpat,ko,'N',Lud,LSigd,wod);
    [ldO, ldOd]     = omega_ar_logdet    (Lu, LSig, jpat, Lud, LSigd);
    if MISS
      [Rlam, Rlamd] = profile_ata        (Laomh, ko, km, p, p, Laomhd);
      [P, Pd]       = profile_mult       (Laomh, Vhat,ko,km,p,p,q,Laomhd,Vhatd);
      if MEAN
        Laomhd        = add_mu_deriv       (r, Laomhd);
      end
      [Lw, Lwd]     = profile_mult       (Laomh, woh, ko, km,p,p,n,Laomhd,wohd);
      clear                              ('Laomh', 'Laomhd')
      [RV, RVd]     = profile_ata        (Vhat, ko, km, p, q, Vhatd);
      if MEAN
        Vhatd         = add_mu_deriv       (r, Vhatd);
      end                                
      [Vw, Vwd]     = profile_mult       (Vhat, woh, ko, km, p, q,n,Vhatd,wohd);
      clear                              ('Vhat', 'Vhatd')
      [Sm, Smd]     = find_Sm            (Scol, miss, Scold);
      [SmP, SmPd]   = atb_deriv          (Sm, Smd, P, Pd);       clear P Pd
      R             = Sm + RV - SmP - SmP';
      Rd            = Smd + RVd - SmPd - permute (SmPd,[2,1,3]);
      [R, Rd]       = atba_c             (Sm, Rlam, R, Smd...
        ,                                Rlamd, Rd);            clear RLam RLamd
      [LR, LRd]     = chol_deriv         (R, Rd);
      [K, Kd]       = forward_sub_deriv  (LR, LRd, RV-SmP...
        ,                                RVd-SmPd);             clear SmP SmPd
      Q             = Sm - RV + K'*K;                           clear RV
      Qd            = ata_deriv          (K, Kd) + Smd - RVd;   clear RVd
      [LQ, LQd]     = chol_deriv         (Q, Qd);               clear Q Qd
      if MEAN
        [Smd, Kd]     = add_mu_deriv       (r, Smd, Kd);
        [LRd, LQd]    = add_mu_deriv       (r, LRd, LQd);
      end
      VwSmLw        = Vw - Sm'*Lw;
      VwSmLwd       = Vwd - atb_deriv    (Sm, Smd, Lw, Lwd);
      [u, ud]       = forward_sub_deriv  (LR, LRd, VwSmLw, VwSmLwd);
      VwKu          = Vw - K'*u;
      VwKud         = Vwd - atb_deriv    (K, Kd, u, ud);
      [v, vd]       = forward_sub_deriv  (LQ, LQd, VwKu, VwKud);
      %
      xSxo = woh'*woh - u'*u + v'*v;
      xSxod = 2*(woh'*Jac(wohd) - u'*Jac(ud) + v'*Jac(vd));
      %
      [ldLR, ldLRd]  = logdet_L          (LR, LRd);
      [ldLQ, ldLQd]  = logdet_L          (LQ, LQd);
      [LSm, LSmd]    = chol_deriv        (Sm, Smd);
      [ldSm, ldSmd]  = logdet_L          (LSm, LSmd);
      ldSo           = ldO + 2*ldLR + 2*ldLQ - 4*ldSm;
      ldSod          = ldOd + 2*ldLRd + 2*ldLQd - 4*ldSmd;
    else
      xSxo  = woh'*woh;
      xSxod = 2*(woh'*Jac(wohd));
      ldSo  = ldO;
      ldSod = ldOd;
    end
    %
    ll = -0.5*(nObs*log(2*pi) + ldSo + xSxo);    
    lld = -0.5*(ldSod + xSxod);
    varargout{1} = lld;
  end
end

function J = Jac(xd) % Change 3-dim derivative to Jacobian matrix
  J = permute(xd, [1,3,2]);
end
