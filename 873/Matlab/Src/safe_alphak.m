%
%   File:              safe_alphak.m
%
%   Description:       Safeguards alphak
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
%   Functions called:  disp, sprintf
%
%   Call: alpha = ...
%         safe_alphak(alpha,alphaL,alphaU,deltaU,ipointa,ipointb,message_level);
%

function [alpha] = ...
          safe_alphak(alpha,alphaL,alphaU,deltaU,ipointa,ipointb,message_level)

    alpha1  = alpha;

    lambdaa = ipointa.lambda; 
    fia     = ipointa.fi;
    norxa   = ipointa.norx;

    lambdab = ipointb.lambda; 
    fib     = ipointb.fi;
    norxb   = ipointb.norx;

    bisection = 0;

    if (alpha < alphaL | alpha > alphaU)
       if (norxb < norxa)
          alpha = deltaU + fib + (norxb^2)*(deltaU - lambdab);
       else
          alpha = deltaU + fia + (norxa^2)*(deltaU - lambdaa);
       end
       if (alpha < alphaL | alpha > alphaU)
          bisection = 1;
          alpha     = (alphaL + alphaU)/2;
       end
    end

    if (message_level == 2)
       disp('-----------')
       disp('safe_alphak')
       disp('-----------')
       if (alpha ~= alpha1)
          disp(sprintf('\nsafeguarding alphak'))
          if (bisection)
             disp(sprintf('using bisection!!'))
          end 
       else
          disp(sprintf('\nnot safeguarding alphak'))
       end
       disp(sprintf('alpha = %e\n',alpha))
    end
