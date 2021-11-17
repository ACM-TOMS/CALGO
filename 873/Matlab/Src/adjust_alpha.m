%
%   File:              adjust_alpha.m
%
%   Description:       Adjusts the parameter alpha until a suitable
%                      interpolation point can be computed
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
%   Functions called:  abs, any, b_epairs, disp, max, smallnu, sprintf
%
%   Call: [alpha,alphaU,nconv,epair1,epair2,z,deltaU,adjk,iterate,v1] =        ...
%          adjust_alpha(nconv,epair1,epair2,Balpha,norg,alphaL,alphaU,deltaU,  ...
%                       epsilon,z,eigensolver,eigensolverpar,lopts);
%

function [alpha,alphaU,nconv,epair1,epair2,z,deltaU,adjk,iterate,v1] =       ...
          adjust_alpha(nconv,epair1,epair2,Balpha,norg,alphaL,alphaU,deltaU, ...
                       epsilon,z,eigensolver,eigensolverpar,lopts)

     v1      = eigensolverpar.v0;

     iterate = 0;

     if (nconv > 0)
        adjust = smallnu(epair1.anu,epair1.noru,norg,epsilon.nu);
        if (~adjust)
           iterate = 1;
        elseif (~any(z))
           z = epair1.u;
        end
        if (nconv > 1 & adjust)
           adjust = smallnu(epair2.anu,epair2.noru,norg,epsilon.nu);
           if (~adjust)
              iterate = 2;
           end
        end
     else
        adjust = 1;
     end

     adjk = 0;

     smallalphaint = abs(alphaU-alphaL) <= ...
                     epsilon.alpha*max(abs(alphaL),abs(alphaU));

     alpha = Balpha.alpha;

     while (adjust & ~smallalphaint)

           alphaU = alpha;
           alpha  = (alphaL + alphaU)/2;

           Balpha.alpha  = alpha;

           smallalphaint = abs(alphaU-alphaL) <= ...
                           epsilon.alpha*max(abs(alphaL),abs(alphaU));

           [nconv,epair1,epair2,v1] = ...
            b_epairs(Balpha,eigensolver,eigensolverpar,lopts);
           eigensolverpar.v0 = v1;

            adjk = adjk + 1;

            if (nconv > 0)
               deltaU = upd_deltaU(epair1,Balpha.g,deltaU,lopts.message_level);
               adjust = smallnu(epair1.anu,epair1.noru,norg,epsilon.nu);
               if (~adjust)
                  iterate = 1;
               elseif (~any(z))
                  z = epair1.u;
               end
               if (nconv > 1 & adjust)
                  adjust  = smallnu(epair2.anu,epair2.noru,norg,epsilon.nu);
                  if (~adjust)
                     iterate = 2;
                  end
               end
            else
               adjust = 1;
            end
     end

     if (lopts.message_level == 2)
        disp('------------')
        disp('adjust_alpha')
        disp('------------')
        disp(sprintf('\n# of times alpha was adjusted: %d\n',adjk))
     end

