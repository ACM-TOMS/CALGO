%
%   File:              interpol2.m
%
%   Description:       2-point interpolation scheme 
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
%   Functions called:  abs, disp, sprintf
%
%   Call: alpha = interpol2(ipointa,ipointb,deltaU,Delta,alphaL,alphaU,message_level);
%

function [alpha] = interpol2(ipointa,ipointb,deltaU,Delta,alphaL,alphaU,message_level)

     if (message_level == 2)
        disp('---------');
        disp('interpol2');
        disp('---------');
        disp(' ');
     end

     lambdaa   = ipointa.lambda;
     fia       = ipointa.fi;
     norxa     = ipointa.norx;

     lambdab   = ipointb.lambda;
     fib       = ipointb.fi;
     norxb     = ipointb.norx;

     den    = Delta * (norxb - norxa);
     mineps = 2e-308;
     if (abs(den) <= mineps)
        lambdabar = deltaU;
        if (message_level == 2)
           disp(sprintf('safeguarding lambdabar, divide by zero\n'))
        end
     else
        lambdabar = (lambdaa * norxa * (norxb - Delta) + ...
                     lambdab * norxb * (Delta - norxa)) / den;
        if (lambdabar > deltaU)
           lambdabar = deltaU;
           if (message_level == 2)
              disp(sprintf('safeguarding lambdabar, greater than deltaU\n'))
           end
        end
     end

     alphabara = lambdaa + fia;
     alphabarb = lambdab + fib;

     if (abs(lambdab - lambdaa) <= mineps)
        alpha = (alphaL + alphaU)/2;
        if (message_level == 2)
           disp(sprintf('Cannot compute omega. Using bisection\n'));
        end
     else
        omega   = (lambdab - lambdabar)/(lambdab - lambdaa);
        omegaco = 1 - omega;

        num     = norxa*norxb*(norxb - norxa)*(lambdaa - lambdabar)* ...
                  (lambdab - lambdabar);
        den     = (omega*norxb + omegaco*norxa)*(lambdab - lambdaa);
        if (abs(den) <= mineps)
           alpha = (alphaL + alphaU)/2;
           if (message_level == 2)
              disp(sprintf('Cannot compute alpha. Using bisection\n'));
           end
        else
           alpha = omega * alphabara + omegaco*alphabarb + num/den; 
        end
     end


