% FREUDPOLYNOMIALS
%
% References: (1)  J. of Scientific Computing, 	v. 30, no. 3, (2007) p. 409-440
%             (2)  Applied and Computational Harmonic Analysis, v. 20, no. 1, (2006) 3-25

% Inputs:
%    x         uniformly spaced points on [-1,1) which serve as trapezoid rule quadrature points
%              the wieghts and Freud polynomials are also evaluated on this grid
%    N         evaluates Freud polynomial orders 0,1,...,N on grid x
%    lambda    exponent parameter in the weight functions
%    small     small number, e.g. 10e-24 <= small <= eps, see references (1) and (2) for details
%    g         (optional)  if the polynomial values are disered at one or very few points, freudPolynomials
%              should be called once with a larger number of x to calculate the recurrence coefficients g
%              and then called again passing g as an input
% Outputs:
%    psi       Freud polynomial orders 0,1,...,N on grid x
%      w       weights
%      g       recurrence coefficients
% Called by:
%    frp.m
% Last modified: October 17, 2007


function [psi,w,g] = freudPolynomials(x,N,lambda,small,g)
 
    M = length(x);
    c = -log(small);
    lambda = floor(lambda);
    w = exp(-c.*x.^(2*lambda));                                    % weights

    psi = zeros(N+1,M);
    psi(1,:) = ones(1,M);
    psi(2,:) = x;
    
    if nargin<5                                                     % compute the recursion coefficients 
                                                                    % and evaluate the polynomials
           g = zeros(1,N+1);
        g(1) = 2*mean(w.*psi(1,:).^2);  
        g(2) = 2*mean(w.*psi(2,:).^2);      
                    
         for k=2:N
           psi(k+1,:) = x.*psi(k,:) - g(k).*psi(k-1,:)./g(k-1);   % equation (2.14) of ref. (1)
               g(k+1) = 2*mean(w.*psi(k+1,:).^2);                 % Trapezoid rule
         end
         
    else                                                          % given the recursion coefficients, evaluate
                                                                  % the polynomials
         for k=2:N
           psi(k+1,:) = x.*psi(k,:) - g(k).*psi(k-1,:)./g(k-1); 
         end
    
    end
  
  
