% CHEBYSHEVPADE  Rational reconstruction of Chebyshev approximations
%
% References: (1) C.K. Clenshaw and K. Lord, Rational approximations from Chebyshev
%                  series, in: Studies in Numerical Analysis, ed. B.K.P. Scaife, Academic
%                  Press, 1974, pp. 95–113.
% Inputs
%    coeff    Chebyshev coefficients
%    M      (M+1) terms in the denominator
%    Nc     only length(coeff) - Nc Chebyshev coefficients are
%              used to build the rational approximation
%    xp     grid to evaluate the rational approximant on
% Outputs
%    up     the rational approximation evaluated on grid xp
% Notes:
%     Linear Chebyshev Pade approximation (see reference (1))
%      rT(x) = (p0*T0 + p1*T1 +...+ pn*Tn) / (q0*T0 + q1*T1 +...+ qm*Tm)
% Called by:
%   1) postProcessDriver.m,  2) chebyshevPade_example.m
% Example usage:
%       see chebyshevPade_example.m
% Last modified: October 17, 2007

  function up = chebyshevPade(coeff,M,Nc,xp)

  coeff(1) = 2*coeff(1);      % undo halving of first and last coeff
  coeff(end) = 2*coeff(end);

    Na = length(coeff);
    Na = Na - Nc;             % number of coefficients to be used

    N = Na - 2*M - 1;
  
    K = N + M;

    P = zeros(1,N+1);
    Q = zeros(1,M+1);
    A = zeros(M,M);
    b = zeros(M,1);

% ------------ find q's -----------------------------------

  ki = 1;
  for k = N+1:K
     b(ki) = -0.5*( coeff(k+0+1)+coeff(abs(k-0)+1) );
     for m = 1:M
       A(ki,m) = 0.5*( coeff(k+m+1)+coeff(abs(k-m)+1) );
     end
     ki = ki + 1;
  end

  Q(1) = 1.0;                    % fix q_0 = 1
  Q(2:M+1) = A\b;                % solve linear system to find q_1,...,q_M

% ------------ find p's -----------------------------------

  for n = 0:N
    P(n+1) = coeff(n+1);
    for m = 1:M
      P(n+1) = P(n+1) + 0.5*( coeff(m+n+1) + coeff(abs(m-n)+1) )*Q(m+1);
    end
    if n==0
      P(n+1) = 0.5*P(n+1);            % halve first p coefficient p0
    end
  end

% ------- evaluate Chebyshev-Pade approximant on xp grid ---------------

  up = zeros(length(xp),1);

  for k=1:length(xp)

    num = 0;
    for j=0:length(P)-1
      num = num + P(j+1)*cos(j*acos(xp(k)));
    end

    dem = 0;
    for j=0:length(Q)-1
      dem = dem + Q(j+1)*cos(j*acos(xp(k)));
    end

    up(k) = num/dem;

  end

% ----------------------------------------------------------------------
