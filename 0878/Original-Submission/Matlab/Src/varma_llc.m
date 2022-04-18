% VARMA_LLC  VARMA likelihood and its derivative for complete data
%
%   [LL, OK] = VARMA_LLC(X, A, B, Sig) sets LL to the exact log-likelihood
%   function value for n complete observations of the r-variate zero-mean VARMA
%   time series:
%
%            x(t) = A1¸x(t-1) + ... + Ap¸x(t-p) + y(t)
%   where
%            y(t) = eps(t) + B1¸eps(t-1) + ... + Bq¸eps(t-q),
%
%   x(t) is r-dimensional and eps(t) is r-variate normal with mean 0 and
%   covariance Sig. Sig and the Ai's and Bi's are rÚr matrices. On entry A
%   should be the r Ú r¸p matrix [A1 A2...Ap], B should be [B1 B2...Bq] and X
%   should have x(t) in its t-th column for t = 1,...,n. OK = true indicates
%   success, but is false if the model is non-stationary.
%
%   [LL, OK, LLD] = VARMA_LLC(X, A, B, Sig) returns the log-likelihood function
%   value in LL and its gradient in LLD.
%
%   [LL, OK, EPS] = VARMA_LLC(X, A, B, Sig, 'res') returns the maximum
%   likelihood estimate of the residuals in EPS.
%
%   [LL, OK, LLD] = VARMA_LLC(X, A, B, Sig, J) may be used to speed up the
%   gradient calculation when A, B and/or Sig depend on a smaller set of
%   independent variables, and only the gradient w.r.t. this smaller set is
%   sought (for instance to fit structural models, distributed lags, and other
%   models with constraints on the parameter matrices). J gives the Jacobian of
%   the change of variables and can have r^2¸p, r^2¸(p+q) or r^2¸(p+q) +
%   r¸(r+1)/2 rows, for when A, A and B, or A, B and Sig depend on a smaller
%   set, respectively (the column count of J equals the number of variables in
%   the smaller set). Thus the possible variable changes are: theta --> A, theta
%   --> [A B], or, theta --> [A B vech(Sig)].
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

function [ll, ok, out3] = varma_llc(X, A, B, Sig, J_code)
  code = ''; J = [];
  if nargin > 4, if ischar(J_code), code = J_code; else J = J_code; end, end
  FINDRES = isequal(code,'res');
  DIFF = nargout == 3 && ~FINDRES;
  CHGVAR = nargin>=5 && ~FINDRES;
  [p, q, r, n] = get_dimensions(A, B, Sig, X);
  ko = 0:r:r*n;
  x = X(:);
  ll = 0; out3 = 0;
  PLU = vyw_factorize(A);
  vyw_ok = isempty(PLU) || isempty(PLU{1}) || PLU{1}(1)~=0;
  ok = vyw_ok && is_stationary(A, Sig, PLU);
  if ~ok, 
    if nargout<=1, error('Non-stationary model'); else return, end
  end
  %
  if nargout <= 2 || FINDRES % ONLY FUNCTION VALUE
    [C, G, W]      = find_CGW          (A, B, Sig);
    S              = vyw_solve         (A, PLU, G);
    [Su, Olow]     = omega_build       (S, G, W, p, n);
    [Lu, Ll, info] = omega_factor      (Su, Olow, p, q, ko); ascertain(info==0);
    w              = lambda_multiply   (A, x, false(r, n));
    z              = omega_forward     (Lu, Ll, w, p, q, ko);
    ld             = omega_logdet      (Lu, Ll, p, q, ko);

    ll = -(r*n*log(2*pi) + z'*z + ld)/2;

    if FINDRES
      eps = res_miss(A, C, Lu, Ll, w);
      out3 = eps;
    end

  else  % FIND ALSO DERIVATIVES
    [CCd,GGd,WWd] = find_CGW_deriv      (A, B, Sig);
    S             = vyw_solve           (A, PLU, GGd);
    RHS           = vyw_deriv_rhs       (A, GGd, S);
    [G, Gd]       = der2array           (GGd);
    [W, Wd]       = der2array           (WWd);
    if CHGVAR
      RHS = chng_var (RHS, J);
      Gd  = chng_var (Gd, J);
      Wd  = chng_var (Wd, J);
    end
    Sd            = vyw_solve           (A, PLU, RHS);
    [Su, Olow]    = omega_build         (S, G, W, p, n);
    [Sud, Olowd]  = omega_build_deriv   (Sd, Gd, Wd, p, n);
    [Lu, Ll,info] = omega_factor        (Su, Olow, p, q, ko); ascertain(info==0);
    [Lud, Lld]    = omega_factor_deriv  (Sud, Olowd, Lu, Ll, p, q, ko);
    xd            = zeros               (r*n, 1, r^2*(p+q)+r*(r+1)/2);
    [w, wd]       = lambda_multiply     (A, x, false(r, n), xd);
    if CHGVAR, wd = chng_var            (wd, J); end
    z             = omega_forward       (Lu, Ll, w, p, q, ko);
    zd            = omega_forward_deriv (Lu, Ll, Lud, Lld, z, wd, p, q, ko);
    [logd,logdd]  = omega_logdet        (Lu, Ll, p, q, ko, Lud, Lld);
    
    ll = -(r*n*log(2*pi) + z'*z + logd)/2;
    lld = -z'*squeeze(zd) - logdd/2;
    out3 = lld;
  end
end
