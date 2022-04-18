function [p,sigma] = ms(g,delta,gamma_qN,tolmin);
%
%  Implementation of the More'-Sorensen method to solve
%
%          min g'p + \half p'Bp   s.t. ||p|| <= delta
%
%  where B is an L-BFGS matrix with B_0 = (1/gamma_qn)*I 
%


global S Y RHO;

itn     = 0;
max_itn = 500;
sigma   = 0;
epsil   = sqrt(eps);   %lower bound on gamma*sigma
show    = 1;

%solve Bp = -g
[n,stored_vecs] = size(S);

p     = two_loop(-g,gamma_qN,stored_vecs);
Bp    = -g;
pnorm = norm(p);

if ~((pnorm<=delta) | (abs(pnorm-delta) <= tolmin*delta ))

  %form B for cholesky decomposition

  if stored_vecs > 0
    [a_unroll,b_unroll]=unrolling(1/gamma_qN);
  end
  
  for itn=1:max_itn
      
    B     =  qN_B(a_unroll,b_unroll,gamma_qN,stored_vecs);
    R     =  chol(B+sigma*eye(n));
    z     = -R'\g;  
    p     =  R\z; 
    q     =  R'\p; 
    pnorm =  norm(p);
    qnorm =  norm(q);

    % update sigma
    sigma = sigma + ((pnorm-delta)*(pnorm)^2)/(delta*(qnorm^2));
    
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

