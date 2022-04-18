%
%   File:              convergence.m
%
%   Description:       Check conditions for convergence
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
%   Date:              15 March 2007
%    
%   Functions called:  abs, boundary_sol, disp, interior_sol, max, norm, 
%                      quasioptimal_sol, sprintf
%
%   Call:  [convsol,xtilde,lambdatilde,stopcond] =                         ...
%           convergence(Balpha,nconv,epair1,epair2,Delta,Deltasqp1,        ...
%                       alphaL,alphaU,iter,maxiter,iterate,epsilon,message_level)
%

function [convsol,xtilde,lambdatilde,stopcond] =                         ...
          convergence(Balpha,nconv,epair1,epair2,Delta,Deltasqp1,        ...
                      alphaL,alphaU,iter,maxiter,iterate,epsilon,message_level)

   if (message_level >= 1)
      disp(sprintf('LSTRS iteration: %d',iter));
   end

   if (iterate == 0)
      boundsol = 0;
   else
      lambda1 = epair1.lambda;
      if (iterate == 1)
         trsiterate = epair1;
      else
         trsiterate = epair2;
         if (message_level == 2)
            disp('Using 2nd eigenvalue')
         end
      end
      boundsol = boundary_sol(trsiterate,lambda1,Delta,epsilon.Delta,message_level);
   end

   if (nconv > 0)
      intsol = interior_sol(epair1,Delta,epsilon.Int,message_level);
   else
      intsol = 0;
   end

   [qopsol,xtilde,lambdatilde] = ...
    quasioptimal_sol(Balpha,nconv,epair1,epair2,Deltasqp1,epsilon.HC,message_level);

   if (iterate == 0)
      smallalphaint = 1;
   else
      smallalphaint = abs(alphaU-alphaL) <= epsilon.alpha*max(abs(alphaL),abs(alphaU));
   end

   convsol = boundsol | intsol | qopsol | smallalphaint | iter > maxiter;

   stopcond.boundsol      = boundsol;
   stopcond.intsol        = intsol;
   stopcond.qopsol        = qopsol;
   stopcond.smallalphaint = smallalphaint;
   stopcond.iterlim       = iter > maxiter;

   if (message_level >= 1)
      if (iterate ~= 0)
         x = trsiterate.u/trsiterate.nu;
         relerr = abs(norm(x)-Delta)/Delta;
         disp(sprintf('||x||: %e, lambda: %e',norm(x),lambda1))
         disp(sprintf('|||x||-Delta|/Delta:     %e\n',relerr));
      end
   end

   if (message_level == 2)
      disp('------------');
      disp('convergence');
      disp('------------');
      disp(sprintf('\niter:  %d',iter));
      disp(sprintf('boundary solution:        %d',boundsol));
      disp(sprintf('interior solution:        %d',intsol));
      disp(sprintf('quasi-optimal solution:   %d',qopsol));
      disp(sprintf('small alpha interval:     %d',smallalphaint));
      disp(sprintf('max number of iterations: %d\n',iter > maxiter));
   end
