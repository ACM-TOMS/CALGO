%
%   File:              adjust_eigentol.m
%
%   Description:       Adjusts the desired accuracy in the eigenpairs
%                      according to the relative accuracy in the solution.
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
%   Functions called:  abs, disp, max, min, norm, sprintf
%
%   Call: [eigentol] = ...
%          adjust_eigentol(eigentol,iterate,epair1,epair2,Delta,maxeigentol,message_level);
%

function [eigentol] = ... 
         adjust_eigentol(eigentol,iterate,epair1,epair2,Delta,maxeigentol,message_level);

     if (iterate ~= 0)
        if (iterate == 1)
           x = epair1.u/epair1.nu;
        elseif (iterate == 2)
           x = epair2.u/epair2.nu;
        end
        relerr = abs(Delta-norm(x))/Delta;
        if (isa(maxeigentol,'double'))
           eigentol = max( min(eigentol,relerr), maxeigentol );
        elseif (relerr < maxeigentol.itermaxacc)
           eigentol = maxeigentol.maxeigentol;
        end
     end

     if (message_level == 2)
        disp('--------------')
        disp('adjust_eig_tol')
        disp('--------------')
        disp(sprintf('\nadjusted eigentol = %6.4e',eigentol))
     end
