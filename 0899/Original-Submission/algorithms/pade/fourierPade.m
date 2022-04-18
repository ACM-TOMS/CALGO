% FOURIERPADE  Rational reconstruction of Fourier approximations 
%
%  References: 1) Fourier-Padé approximations and filtering for spectral simulations of an incompressible Boussinesq
%                 convection problem, Mathematics of Computation, v. 76, 2007, p. 1275-1290.
%
% Inputs
%    uh    Fourier coefficients
%    M    (M+1) terms in the denominator
%    Nc    only length(coeff) - Nc Fourier coefficients are
%              used to build the rational approximation
%    xp    grid to evaluate the rational approximant on
% Outputs
%     R    rational approximation R = (P_K)/(Q_M)
%  Called by:
%          1) postProcessDriver.m
% Last modified: October 17, 2007

 function R = fourierPade(uh,M,Nc,xp)

     a =  -1;  b =  1;
     xp = pi + 2*pi*(xp-b)/(b-a);     % map [a,b] to [-pi,pi]

      N2 = length(uh);
      uh = ((-1).^(0:N2-1)).*fftshift(uh)/N2;    % reorder Foruier coefficients
      N = round(N2/2);

     j = 0:N2-1;
     x = -pi + j*pi./(N);
     dx = 1.0/(2*N);
     K = N - 2*M;
     K = K - Nc;
     Nh = N + 1;

% --- form matrix A, (3.9), and -----------------------------------

     ni = 0;
     for n=-K-M:-K-1
        ni = ni + 1;
        mi = 0;
        for m=-M:M
          mi = mi + 1;
          A(ni,mi) = uh(Nh+n-m);
          cm(ni,mi) = n-m;
        end
     end

     for n=K:K-1+M
        ni = ni + 1;
        mi = 0;
        for m=-M:M
          mi = mi + 1;
          A(ni,mi) = uh(Nh+n-m);
          cm(ni,mi) = n-m;
        end
     end

% -- solve the linear system Ac = 0 to find the coefficients of the polynomial
% -- Q_M in the demoninator

     c = null(A);

% -- (3.6), calculate the coefficients of the numerator, P_K

     kAdj = K+1;
     mAdj = M+1;
     for k=-K:K-1
       b(k+kAdj) = 0;
       for m=-M:M
         b(k+kAdj) = b(k+kAdj) + c(m+mAdj)*uh(Nh+k-m);
       end
     end

% -- (3.5), form P and Q

     for j=1:length(xp)
       P(j) = 0;
       for k=-K:K-1
         P(j) = P(j) + b(k+kAdj)*exp(i*k*xp(j));
       end
     end

     for j=1:length(xp)
       Q(j) = 0;
       for m=-M:M
         Q(j) = Q(j) + c(m+mAdj)*exp(i*m*xp(j));
       end
     end

   R = real( P./Q );
