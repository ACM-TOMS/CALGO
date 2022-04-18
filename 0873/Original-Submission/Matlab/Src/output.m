%
%   File:              output.m
%
%   Description:       Sets output parameters according to the
%                      stopping criteria that were satisfied
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
%   Date:              13 April 2007
%
%   Functions called:  any, cg, correct, disp, func2str, isa, length, matvec, 
%                      norm, smallnu, sprintf, strcmp
%
%   Call: [x,lambda,info,moreinfo] =                                     ...
%          output(epair1,epair2,iterate,lambdatilde,xtilde,Delta,        ...
%                 Balpha,norg,z,iter,adj,stopcond,lopts);
%

function [x,lambda,info,moreinfo] =                                      ...
          output(epair1,epair2,iterate,lambdatilde,xtilde,Delta,         ...
                 Balpha,norg,z,iter,adj,stopcond,lopts)

%
% Declare global variable to count matrix-vector products
% Routines that use it: lstrs_method, matvec, output
%
global mvp_lstrs;

   boundsol      = stopcond.boundsol;
   intsol        = stopcond.intsol;
   qopsol        = stopcond.qopsol;
   smallalphaint = stopcond.smallalphaint;
   iterlim       = stopcond.iterlim;

   clear stopcond;

   nstopcond = 5;
   for i = 1:nstopcond 
      stopcond(i).on = 0;
   end
   stopcond(1).message = 'Boundary Solution';
   stopcond(1).code    = 'BS';
   stopcond(2).message = 'Interior Solution';
   stopcond(2).code    = 'IS';
   stopcond(3).message = 'Quasi-optimal Solution';
   stopcond(3).code    = 'QO';
   stopcond(4).message = 'Safeguarding interval for alpha is too small';
   stopcond(4).code    = 'AI';
   stopcond(5).message = 'Maximum number of iterations reached';
   stopcond(5).code    = 'IT';
%   stopcond(6).message = 'KKT condition was satisfied';


   message_level = lopts.message_level;
  
   havesol   = 0;

   Balpha.bord = 0;
   n           = Balpha.dim - 1;

   if (message_level == 2)
      disp('-------');
      disp('output');
      disp('-------');
      if (~isempty(epair1))
         disp(sprintf('lambda1: %e',epair1.lambda));
         disp(sprintf('nu1:     %e',epair1.nu));
         if (~isempty(epair2))
            disp(sprintf('lambda2: %e',epair2.lambda));
            disp(sprintf('nu2:     %e\n',epair2.nu));
         end
      end
   end

   if (boundsol)
      stopcond(1).on = 1;
      havesol        = 1;
      info           = 0;
      lambda         = epair1.lambda;
      if (iterate == 1)
         x = epair1.u/epair1.nu;
      else
         x = epair2.u/epair2.nu;
         if (message_level == 2)
            disp('Terminates with SECOND eigenvalue!');
         end
      end
   end

   if (intsol)
      stopcond(2).on = 1;
      if (~havesol)
         havesol = 1;
         if (strcmp(lower(lopts.interior(1)),'y'))
            if (message_level > 0)
               disp(sprintf('Computing an Interior Solution\n'));
               disp(sprintf('Solving the system: Hx = -g with CG ...\n'));
            end

            info   = 1;

            x      = cg(Balpha,-Balpha.g,lopts.intsoltol);
            lambda = 0;
         else
            lambda = epair1.lambda;
            info   = -1;
            if (iterate == 1)
               x = epair1.u/epair1.nu;
            else
               x = epair2.u/epair2.nu;
               if (message_level == 2)
                  disp('Terminates with SECOND eigenvalue!');
               end
            end
         end
      end
   end

   if (qopsol)
      stopcond(3).on = 1;
      if (~havesol)
         havesol        = 1;
         info           = 2;
         lambda         = lambdatilde;
         x              = xtilde;
      end
   end

   if (smallalphaint)
      stopcond(4).on = 1;
      if (~havesol)
         havesol = 1;
         lambda  = epair1.lambda;

         if (iterate == 0)
            info = -4;
            x    = [];
         else   
            info = -2;
            if (iterate == 1)
               x = epair1.u/epair1.nu;
            else
               x = epair2.u/epair2.nu;
            end
            norx = norm(x);
            if (norx < Delta)
               if (message_level == 2)
                  disp('Hard Case !!')
               end
               if (any(z) & strcmp(lower(lopts.correction(1)),'y'))
                  x = correct(z,x,norx,Delta);
                  if (message_level == 2)
                     disp('A correction term was added to x.')
                  end
               end
            end
         end
      end
   end

   if (iterlim)
      stopcond(5).on = 1;
      if (~havesol)
         lambda = epair1.lambda;
         if (iterate == 0)
            info = -4;
            x    = [];
         else   
            info = -3;
            if (iterate == 1)
               x = epair1.u/epair1.nu;
            else
               x = epair2.u/epair2.nu;
               if (message_level == 2)
                  disp('Terminates with SECOND eigenvalue!');
               end
            end
         end
      end
   end

   if (iterate == 0)
      kkt = -1;
   else
      kkt = Balpha.g + matvec(x,Balpha) - lambda*x;
      kkt = norm(kkt)/norg;
   end

   moreinfo.iter   = iter;
   moreinfo.solves = adj;
   moreinfo.kkt    = kkt;
   moreinfo.alpha  = Balpha.alpha;
%
% Record value of global counter of number of matrix-vector products
%
   moreinfo.mvp    = mvp_lstrs; 

   if (message_level > 0)
      disp(sprintf('\nNumber of LSTRS Iterations:     %d',moreinfo.iter))
      disp(sprintf('Number of calls to eigensolver: %d',moreinfo.solves))
      disp(sprintf('Number of MV products:          %d\n',moreinfo.mvp))

      if (iterate ~= 0)
         disp(sprintf('\n(||x||-Delta)/Delta: %e\n',abs(norm(x)-Delta)/Delta));
         disp(sprintf('lambda: %e\n',lambda));
         disp(sprintf('||g + (H-lambda* I)x||/||g|| = %e\n',kkt));
      end

      i = 1;
      while (~stopcond(i).on)
            i = i+1;
      end
      exitcond = stopcond(i).code;
      if (i <= 3)
         if (i ~= 2)
            disp(sprintf('The vector x is a %s',stopcond(i).message));
         else
            if (strcmp(lower(lopts.interior(1)),'y'))
               disp(sprintf('The vector x is an %s',stopcond(i).message));
            else
               disp(sprintf('The vector x is the last iterate of LSTRS'));
            end
         end
      else
         if (info ~= -4)
            disp(sprintf('The vector x is the last iterate of LSTRS'))
         else
            disp(sprintf('The vector x is empty: no iterate is available'))
         end
         disp(sprintf('LSTRS exit because: %s',stopcond(i).message))
      end

      j = i+1;
      found = 0;
      while (j <= nstopcond & ~stopcond(j).on)
         j = j + 1;
      end             
      found = j <= nstopcond;
         
      if (found)
         disp(sprintf('\nOther Stopping Criteria Satisfied:'))
         i = i+1;
         while (i <= nstopcond)
            if (stopcond(i).on)
               exitcond = [exitcond,',',stopcond(i).code];
               disp(sprintf('   %s',stopcond(i).message))
            end
            i = i+1;
         end
     end
     disp(sprintf('\n'))
  else
      i = 1;
      while (~stopcond(i).on)
            i = i+1;
      end
      exitcond = stopcond(i).code;
      i = i+1;
      while (i <= nstopcond)
         if (stopcond(i).on)
            exitcond = [exitcond,',',stopcond(i).code];
         end
         i = i+1;
      end
   end

   moreinfo.exitcond = exitcond;

   if (info < 0)
      if (info == -4)
         disp('x is empty: no iterate is available!');
      elseif (info == -1)
         disp('An interior solution was detected, but not computed as instructed.');
         disp('If this was a regularization problem, you might want to decrease Delta.');
         disp('x and lambda contain the last iterate available.')
      else
         disp('Could not compute a solution with the desired accuracy.')
         disp('You might want to try the following:')
         disp('- Relax tolerances such as: epsilon.Delta, epsilon.HC, accuracy of eigenpairs')
         disp('- Increase maximum number of iterations allowed')
         disp('x and lambda contain the last iterate available.')
      end
   end


%
% Plots results
%

if (info >= 0 & strcmp(lower(lopts.plot(1)),'y'))
   figure(1);
   clf;
   plot(x,'b','linewidth',3);
   title(sprintf('Problem: %s. Dimension: %d. Delta: %e',lopts.name,n,Delta));
end

