% VARMA_LLM VARMA likelihood for missing data
%
%   [LL, OK] = VARMA_LLM(X, A, B, Sig, mu, miss) returns the value of the exact
%   log-likelihood function for an r-variate VARMA time series of length n:
%
%             x(t) - mu = A1¸(x(t-1) - mu) + ... + Ap¸(x(t-p) - mu) + y(t)
%    where
%             y(t) = eps(t) + B1¸eps(t-1) + ... + Bq¸eps(t-q),
%
%   and x(t), y(t) and eps(t) are r-dimensional with eps(t) normally
%   distributed, N(0,Sig). Sig and the Ai's and Bi's are rÚr matrices, and mu is
%   an r-vector with the mean of x(t). A should be the  r Ú r¸p matrix [A1
%   A2...Ap], B should be [B1 B2...Bq] and X should have x(t) in its t-th column
%   for t = 1,...,n. Set miss(i,t) to true if X(i,t) is missing. OK = true
%   indicates success, but is false if the model is non-stationary.
%
%   [LL, OK, LLD] = VARMA_LLM(X, A, B, Sig, mu, miss) returns the log-likelihood
%   function value in LL and its gradient in LLD.
%
%   [LL, OK, EPS] = VARMA_LLM(...,'res') returns the maximum likelihood estimate
%   of the residuals in EPS and [LL, OK, EPS, XM] = VARMA_LLM(..., 'res_miss')
%   finds also the maximum likelihood estimate of the missing values in XM.
%
%   [LL, OK, LLD] = VARMA_LLM(X, A, B, Sig, mu, miss, J) may be used to speed up
%   the gradient calculation when A, B and/or Sig depend on a smaller set of
%   independent variables, and only the gradient w.r.t. this smaller set is
%   sought (for instance to fit structural models, distributed lags, and other
%   models with constraints on the parameter matrices). J gives the Jacobian of
%   the change of variables and can have r^2¸p, r^2¸(p+q) or r^2¸(p+q) +
%   r¸(r+1)/2 rows, for when A, A and B or A, B and Sig depend on a smaller set,
%   respectively (the column count of J equals the number of variables in the
%   smaller set). Thus the possible variable changes are theta --> A, theta -->
%   [A B] and theta --> [A B vech(Sig)]. The only possibilities for mu are to
%   leave it free or to fix it to be zero. Letting mu = [] results in LL
%   containing the likelihood of a zero mean model and in LLD excluding
%   derivatives w.r.t. mu.
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

function [ll, ok, varargout] = varma_llm(X, A, B, Sig, mu, miss, J_code)
  code = ''; J = [];
  if nargin>6, if ischar(J_code), code = J_code; else J = J_code; end, end
  FINDRES = any(strmatch(code, {'res','res_miss'}, 'exact'));
  CHNGVAR = ~isempty(J);
  mu = mu(:);
  MEAN = ~isempty(mu);
  LowTri = struct('LT', true);
  [p, q, r, n] = get_dimensions(A, B, Sig, X);
  nPar = r^2*(p+q) + r*(r+1)/2;
  [ko,ro,km] = find_missing_info(miss);
  if MEAN, mubar = repmat(mu, n, 1); else mubar = zeros(n*r, 1); end
  x = X(:);
  xo = x(~miss) - mubar(~miss);
  nObs = length(xo);
  ll = 0; varargout = {0,0};
  ll = 0; varargout{1} = 0;
  PLU = vyw_factorize(A);
  vyw_ok = isempty(PLU) || isempty(PLU{1}) || PLU{1}(1)~=0;
  ok = vyw_ok && is_stationary(A, Sig, PLU);
  if ~ok, if nargout<=1, error('Non-stationary model'); else return, end, end
  %
  if nargout<=2 || FINDRES  % GRADIENT NOT REQUESTED
    [C, G, W]        = find_CGW            (A, B, Sig);
    S                = vyw_solve           (A, PLU, G);
    Gneg             = find_Gneg           (A, B, C, n);
    Scol             = S_extend            (A, G, S, n);
    Acol             = find_Acol           (A, r);
    [Su, Olow]       = omega_build         (S, G, W, p, n);
    [Suo, Olowo]     = omega_remove_miss   (Su, Olow, miss);
    [Luo, Llo, info] = omega_ltl           (Suo,Olowo,p,q,ko); ascertain(info==0);
    Laom             = find_lambda_om      (Acol, miss);
    Laomh            = profile_back_sub    (Luo, Llo, Laom, ko, km, p, q, p);
    clear Laom
    V                = find_V              (G, Gneg, Scol, miss);
    Vhat             = profile_back_sub    (Luo, Llo, V, ko, km, p, q, q);
    clear V
    wo               = lambda_multiply     (A, xo, miss);
    wohat            = omega_back_sub      (Luo, Llo, wo, p, q, ko);
    RLam             = profile_ata         (Laomh, ko, km, p, p);
    P                = profile_mult        (Laomh, Vhat, ko, km, p, p, q);
    Lw               = profile_mult        (Laomh, wohat, ko, km, p, p, n);
    if ~FINDRES
      clear Laomh
    end
    RV               = profile_ata         (Vhat, ko, km, p, q);
    Vw               = profile_mult        (Vhat, wohat, ko, km, p, q, n);
    if ~FINDRES
      clear Vhat
    end
    Sm               = find_Sm             (Scol, miss);
    SmP              = Sm*P;                                          clear P
    R                = -SmP-SmP' + atba_c  (Sm, RLam, Sm + RV);       clear Rlam
    LR               = cholf               (R);                       clear R
    K                = linsolve            (LR, RV - SmP, LowTri);    clear SmP
    Q                = Sm - RV + K'*K;                                clear RV
    LQ               = cholf               (Q);                       clear Q
    u                = linsolve            (LR, Vw - Sm*Lw, LowTri);
    v                = linsolve            (LQ, Vw - K'*u, LowTri);
    %
    xSxo = wohat'*wohat - u'*u + v'*v;
    %
    logdetO = omega_logdet(Luo, Llo, p, q, ko);
    logdetLR = sum(log(diag(LR)));
    logdetLQ = sum(log(diag(LQ)));
    logdetSm = 2*sum(log(diag(cholf(Sm))));
    logdetSo = logdetO + 2*logdetLR + 2*logdetLQ - 2*logdetSm;
    %
    ll = -0.5*(nObs*log(2*pi) + logdetSo + xSxo);
    if FINDRES
      [eps,xm] = res_miss(A,C,Luo,Llo,wohat,LR,LQ,K,u,v,Laomh,Vhat,Sm,miss);
      varargout{1} = eps;
      if nargout>=4, varargout{2} = xm + mubar(miss); end
    end
    %
  else  % FIND ALSO GRADIENT
    if MEAN, xod        = xmu_deriv         (nPar + r, miss);
    else     xod        = zeros             (nObs, 1, nPar); end
    [wo, wod]           = lambda_multiply   (A, xo, miss, xod);
    [CCd, GGd, WWd]     = find_CGW_deriv    (A, B, Sig);
    [C, Cd]             = der2array         (CCd);
    [G, Gd]             = der2array         (GGd);
    [W, Wd]             = der2array         (WWd);
    S                   = vyw_solve         (A, PLU, GGd);
    RHS                 = vyw_deriv_rhs     (A, GGd, S);
    Sd                  = vyw_solve         (A, PLU, RHS);
    [Gneg, Gnegd]       = find_Gneg         (A, B, C, n, Cd);
    [Scol, Scold]       = S_extend          (A, G, S, n, Gd, Sd);
    [Acol, Acold]       = find_Acol         (A, r, nPar);
    if CHNGVAR
      wod    = chng_var(wod, J);
      Gd     = chng_var(Gd, J);
      Wd     = chng_var(Wd, J);
      Sd     = chng_var(Sd, J);
      Gnegd  = chng_var(Gnegd, J);
      Scold  = chng_var(Scold, J);
      Acold  = chng_var(Acold, J);
    end
    [Su, Olow]          = omega_build       (S, G, W, p, n);
    [Sud, Olowd]        = omega_build_deriv (Sd, Gd, Wd, p, n);
    [Su,Olow,Sud,Olowd] = omega_remove_miss (Su, Olow, miss, Sud, Olowd);
    [Lu,Ll,flg,Lud,Lld] = omega_ltl         (Su, Olow, p, q, ko, Sud, Olowd);
    clear Su Sud Olow Olowd
    [Laom, Laomd]     = find_lambda_om    (Acol, miss, Acold); % Lambda_om
    [Laomh, Laomhd]   = profile_back_sub  (Lu, Ll, Laom, ko, km, p, q, p...
      ,                                   Lud, Lld, Laomd);
    clear Laom Laomd
    [V, Vd]           = find_V            (G, Gneg, Scol, miss, Gd,Gnegd,Scold);
    [Vhat, Vhatd]     = profile_back_sub  (Lu, Ll, V, ko, km, p, q, q...
      ,                                   Lud, Lld, Vd);
    clear V Vd
    [RLam, RLamd]     = profile_ata       (Laomh, ko, km, p, p, Laomhd);
    [P, Pd]           = profile_mult      (Laomh, Vhat, ko, km, p, p, q...
      ,                                   Laomhd, Vhatd);
    if MEAN
      Laomhd            = add_mu_deriv      (r, Laomhd);
      [Lud, Lld]        = add_mu_deriv      (r, Lud, Lld);
    end
    %
    [woh, wohd]       = omega_back_sub    (Lu, Ll, wo, p, q, ko, Lud, Lld, wod);
    [Lw, Lwd]         = profile_mult      (Laomh, woh, ko,km,p,p,n,Laomhd,wohd);
    clear Laomh Laomhd
    [RV,RVd]          = profile_ata       (Vhat, ko, km, p, q, Vhatd);
    if MEAN
      Vhatd             = add_mu_deriv      (r, Vhatd); 
    end
    [Vw, Vwd]         = profile_mult      (Vhat, woh, ko, km, p,q,n,Vhatd,wohd);
    clear Vhat Vhatd
    [Sm, Smd]         = find_Sm           (Scol, miss, Scold);
    [SmP, SmPd]       = atb_deriv         (Sm, Smd, P, Pd);     clear P Pd
    R                 = Sm + RV - SmP - SmP';
    Rd                = Smd + RVd - SmPd - permute(SmPd,[2,1,3]);
    [R, Rd]           = atba_c            (Sm, RLam, R, Smd...  
      ,                                   RLamd, Rd);           clear RLam RLamd
    [LR, LRd]         = chol_deriv        (R, Rd);              clear R Rd
    [K, Kd]           = forward_sub_deriv (LR, LRd, RV-SmP...
      ,                                   RVd-SmPd);            clear SmP SmPd
    Q                 = Sm - RV + K'*K;                         clear RV
    Qd                = ata_deriv         (K, Kd) + Smd - RVd;  clear RVd
    [LQ, LQd]         = chol_deriv        (Q, Qd);              clear Q Qd
    if MEAN
      [Smd,LRd,LQd,Kd]  = add_mu_deriv      (r, Smd, LRd, LQd, Kd);
    end
    VwSmLw            = Vw - Sm'*Lw;
    VwSmLwd           = Vwd - atb_deriv   (Sm, Smd, Lw, Lwd);
    [u, ud]           = forward_sub_deriv (LR, LRd, VwSmLw, VwSmLwd);
    VwKu              = Vw - K'*u;
    VwKud             = Vwd - atb_deriv   (K, Kd, u, ud);
    [v, vd]           = forward_sub_deriv (LQ, LQd, VwKu, VwKud);
    %
    xSxo = woh'*woh - u'*u + v'*v;
    xSxod = 2*(woh'*Jac(wohd) - u'*Jac(ud) + v'*Jac(vd));
    %
    [ldO, ldOd] = omega_logdet(Lu, Ll, p, q, ko, Lud, Lld);
    [ldLR, ldLRd] = logdet_L(LR, LRd);
    [ldLQ, ldLQd] = logdet_L(LQ, LQd);
    [LSm, LSmd] = chol_deriv(Sm, Smd);
    [ldSm, ldSmd] = logdet_L(LSm, LSmd);
    ldSo = ldO + 2*ldLR + 2*ldLQ - 4*ldSm;
    ldSod = ldOd + 2*ldLRd + 2*ldLQd - 4*ldSmd;
    %
    ll = -0.5*(nObs*log(2*pi) + ldSo + xSxo);    
    lld = -0.5*(ldSod + xSxod);
    varargout{1} = lld;
  end
end

function J = Jac(xd) % Change 3-dim derivative to Jacobian matrix
  J = permute(xd, [1,3,2]);
end
