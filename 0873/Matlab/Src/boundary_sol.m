%
%   File:              boundary_sol.m
%
%   Description:       Test Condition for a Boundary Solution
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
%   Call: boundsol =  ...
%         boundary_sol(trsiterate,lambda1alphak,Delta,eps_Delta,message_level);
%

function [boundsol] = ...
          boundary_sol(trsiterate,lambda1alphak,Delta,eps_Delta,message_level)

     norx      = trsiterate.noru/trsiterate.anu;

     boundsol  = (abs(norx - Delta) <= eps_Delta * Delta) & (lambda1alphak <= 0);

     if (message_level == 2)
        disp('------------')
        disp('boundary_sol')
        disp('------------')
        disp(sprintf('\n  | ||x|| - Delta|     eps_Delta * Delta'))
        disp(sprintf('        %6.4f              %6.4f\n',abs(norx-Delta), ...
                      eps_Delta*Delta))
        disp(sprintf('||x|| = %6.4f  Delta = %6.4f\n',norx,Delta))
        disp(sprintf('lambda1(alphak) = %6.4f\n',lambda1alphak))
     end

