%
%   File:              interior_sol.m
%
%   Description:       Test Condition for an Interior Solution
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
%   Call: intsol = interior_sol(epair,Delta,eps_Int,message_level);
%

function [intsol] = interior_sol(epair,Delta,eps_Int,message_level)

     lambda = epair.lambda;
     absnu  = epair.anu;
     noru   = epair.noru;

     intsol =   (noru < Delta * absnu) & (lambda > - eps_Int);

     if (message_level == 2)
        disp('------------')
        disp('interior_sol')
        disp('------------')
        disp(sprintf('\n    ||u||       Delta*|nu|'))
        disp(sprintf('%e   %6.4f\n',noru,Delta*absnu));
        disp(sprintf('     lambda       -eps_Int'))
        disp(sprintf('%e %e\n',lambda,-eps_Int))
        disp(sprintf('lambda > -eps_Int: %d',lambda > - eps_Int));
        disp(sprintf('||x|| < Delta:     %d\n',noru < Delta * absnu));
     end

