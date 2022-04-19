%
%   File:              quasioptimal_sol.m
%
%   Description:       Test Condition for a quasi-optimal solution
%                      obtained in the Hard Case or nearly Hard Case
%
%   Authors:           Marielba Rojas
%                        mr@imm.dtu.dk
%                      Sandra Santos
%                      Danny Sorensen
%
%   Version:           1.2
%
%   System:            MATLAB 6.0 or higher 
%
%   Date:              3 April 2007
%
%   Functions called:  disp, quadratic, sprintf, sqrt
%
%   Call: [qopsol,xtilde,lambdatilde] = ...
%         quasioptimal_sol(Balpha,nconv,epair1,epair2,Deltasqp1,eps_HC,message_level);
%

function [qopsol,xtilde,lambdatilde] = ...
          quasioptimal_sol(Balpha,nconv,epair1,epair2,Deltasqp1,eps_HC,message_level)


    if (message_level == 2)
       disp('----------------');
       disp('quasioptimal_sol');
       disp('----------------');
       disp(' ');
    end

    eta  = eps_HC/(1 - eps_HC);

    if (nconv >= 2)
       lambda1 = epair1.lambda; 
       nu1     = epair1.nu; 
       lambda2 = epair2.lambda; 
       nu2     = epair2.nu; 

       normfsq   = nu1^2 + nu2^2; 
       normf     = sqrt(normfsq);
       product   = normfsq * Deltasqp1;

       if (product >= 1)
          sqrproductm1 = sqrt(product - 1);
          den          = normfsq * sqrt(Deltasqp1);
          if (product > 1)
             tau1   = (nu1 - nu2 * sqrproductm1)/den;
             tau2   = (nu2 + nu1 * sqrproductm1)/den;
          else
             tau1   = nu1/normf;
             tau2   = nu2/normf;
          end
          u1          = epair1.u;
          u2          = epair2.u;
          xtilde      = (tau1*u1 + tau2*u2)/(tau1*nu1 + tau2*nu2);

          lambdatilde = tau1^2*lambda1 + tau2^2*lambda2;

          quadxtilde  = quadratic(Balpha,xtilde);
          lhsexp1     = (lambda2 - lambda1)*(tau2^2)*Deltasqp1;
          rhsexp1     =  -2*eta*quadxtilde;
          qopsol      =  lhsexp1 <= rhsexp1;

          if (message_level == 2)
             disp(sprintf('tau1 = %6.4e',tau1));
             disp(sprintf('tau2 = %6.4e\n',tau2));
          end
          if (~qopsol & product > 1)
             tau1        = (nu1 + nu2 * sqrproductm1)/den;
             tau2        = (nu2 - nu1 * sqrproductm1)/den;
             xtilde      = (tau1*u1 + tau2*u2)/(tau1*nu1 + tau2*nu2);
             lambdatilde = tau1^2*lambda1 + tau2^2*lambda2;
             quadxtilde  = quadratic(Balpha,xtilde);
             rhsexp2     = -2*eta*quadxtilde;
             lhsexp2     = (lambda2 - lambda1)*(tau2^2)*Deltasqp1;
             qopsol      =  lhsexp2 <= rhsexp2;
             if (message_level == 2)
                disp(sprintf('tau1 = %6.4e',tau1));
                disp(sprintf('tau2 = %6.4e\n',tau2));
             end
          else
             lhsexp2 = 1;
             rhsexp2 = 0;
          end
      else
          qopsol      = 0;
          xtilde      = 0;
          lambdatilde = 0;
          lhsexp1     = 1;
          lhsexp2     = 1;
          rhsexp1     = 0;
          rhsexp2     = 0;
      end
   else
      qopsol      = 0;
      xtilde      = 0;
      lambdatilde = 0;
      lhsexp1     = 1;
      lhsexp2     = 1;
      rhsexp1     = 0;
      rhsexp2     = 0;
   end

   if (message_level == 2)
      if (nconv >= 2)
         disp(sprintf('nu1^2 + nu2^2: %6.4e 1/(1+Delta^2): %6.4e\n', ...
                      normfsq,1/Deltasqp1));
         disp(sprintf('product >= 1: %d\n',product >= 1));
      end
      disp(sprintf('qopsol:       %d\n',qopsol));
      s = '(l2-l1)*tau2^2*(Delta^2+1)';
      disp(sprintf('1. %s: %6.4e, -2*eta*psi(xtilde): %6.4e',s,lhsexp1,rhsexp1));
      disp(sprintf('2. %s: %6.4e, -2*eta*psi(xtilde): %6.4e',s,lhsexp2,rhsexp2));
   end
