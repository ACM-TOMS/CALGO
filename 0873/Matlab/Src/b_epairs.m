%
%   File:              b_epairs.m
%
%   Description:       Computes smallest eigenpairs of bordered matrix
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
%   Functions called:  abs, disp, feval, sprintf, sqrt
%
%   Call: [nconv,epair1,epair2,v1] =      ...
%          b_epairs(Balpha,eigensolver,eigensolverpar,lopts);
%

function [nconv,epair1,epair2,v1] =       ...
          b_epairs(Balpha,eigensolver,eigensolverpar,lopts)

      [nconv,lambda1,y1,lambda2,y2,v1] =  ...
       feval(eigensolver,Balpha,eigensolverpar);
   
      if (lopts.heuristics)
         nconv = 2;
      end

      if (nconv > 0)
         m             = Balpha.dim;
         epair1.lambda = lambda1;
         epair1.nu     = y1(1);
         epair1.anu    = abs(epair1.nu);
         epair1.noru   = sqrt(1 - epair1.nu^2);
         epair1.u      = y1(2:m);
         if (nconv > 1)
            epair2.lambda = lambda2;
            epair2.nu     = y2(1);
            epair2.anu    = abs(epair2.nu);
            epair2.noru   = sqrt(1 - epair2.nu^2);
            epair2.u      = y2(2:m);
         else
            epair2 = [];
         end
      else
         epair1 = [];
         epair2 = [];
      end

     if (lopts.message_level == 2)
        disp('--------')
        disp('b_epairs')
        disp('--------')
        disp(sprintf('\nalpha: %e\n',Balpha.alpha))
        disp(sprintf('# converged eigenvalues: %d',nconv))
        if (nconv > 0)
           disp(sprintf('lambda1: %e',epair1.lambda))
           disp(sprintf('nu1:     %e',epair1.nu))
           if (nconv > 1)
               disp(sprintf('lambda2: %e',epair2.lambda))
               disp(sprintf('nu2:     %e\n',epair2.nu))
           end
        end
     end
