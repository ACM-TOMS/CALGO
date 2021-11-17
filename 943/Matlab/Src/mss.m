function [p,sigma] = mss(g,delta,gamma_qN,tolmin);
% 
%  Implementation of the MSS method to solve
%
%          min g'p + \half p'Bp   s.t. ||p|| <= delta
%
%  where B is an L-BFGS matrix with B_0 = (1/gamma_qn)*I 
%
%  More information on the MSS method can be found in:
%
%     "MSS: MATLAB Software for L-BFGS Trust-Region 
%       Subproblems for Large-Scale Optimization"
%
%  Authors: Jennifer B. Erway and Roummel F. Marcia
%
%  Paper available at 
%
%  http://www.wfu.edu/~erwayjb/publications.html
%
%  Copyright (2012): Jennifer B. Erway and Roummel F. Marcia

global S Y RHO;

itn     = 0;
max_itn = 500;
sigma   = 0;
epsil   = sqrt(eps);   %lower bound on gamma*sigma
show    = 1;

if gamma_qN < epsil
    gamma_qN = sqrt(epsil);
end

%solve Bp = -g
stored_vecs = size(S,2);

p     = two_loop(-g,gamma_qN,stored_vecs);
Bp    = -g;
pnorm = norm(p);

if  ~((pnorm<=delta) | (abs(pnorm-delta) <= tolmin*delta ))
  if stored_vecs > 0
    [a_unroll,b_unroll]=unrolling(1/gamma_qN);
  end
  
  for itn=1:max_itn
      
    % update sigma
    phi = 1/pnorm - 1/delta;
    if sigma<epsil
      gradp = two_loop(-p,gamma_qN,stored_vecs);
      sigma = 0;
    elseif (stored_vecs>0)
      gradp = recursionSolve_qN(a_unroll,b_unroll,1/gamma_qN,sigma,-p);
    else 
      gradp = -(1/(1+sigma))*p;
    end
    
    phi_prime = -(p'*gradp)/((pnorm)^3);
    sigma     = sigma - phi/phi_prime;
    
    if sigma<epsil
      p = two_loop(-g,gamma_qN,stored_vecs);
      sigma = 0;
    elseif stored_vecs>0
      p = recursionSolve_qN(a_unroll,b_unroll,1/gamma_qN,sigma,-g);
    else
      p = -(1/(1+sigma))*g;
    end
    pnorm = norm(p);

    if (abs(pnorm-delta) <= tolmin*delta )
      break;
    end
  end %for loop
end %if loop
  

if show,
  Bp = qN_multiply(a_unroll,b_unroll,p,gamma_qN,stored_vecs);
  error1 = norm(Bp+sigma*p+g);
  error2 = abs(sigma*(sqrt(p'*p)-delta));
  error3 = norm(two_loop(Bp,gamma_qN,stored_vecs)-p); 
  fprintf('Error in 1st optimality condition: %8.2e\n', error1);
  fprintf('Error in 2nd optimality condition: %8.2e\n', error2);
end
