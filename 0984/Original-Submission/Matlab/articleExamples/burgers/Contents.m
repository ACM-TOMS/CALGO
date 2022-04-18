% Burgers ODE Problem
%
% Files
% burgersfun.m           - Burgers ODE function file, taken verbatim from
%                          burgersode.m
% burgersfun_noloop.m    - Burgers ODE function file with loops removed
% burgersRatios_loop.m   - Tests ADiGator, ADiMat, INTLAB, and MAD to 
%                          compute Jacobian with the original ODE function
% burgersRatios_noloop.m - Tests ADiGator, ADiMat, INTLAB, and MAD to 
%                          compute Jacobian with the modified ODE function
% burgerssolve.m         - Solves Burgers ODE using ADiGator and ode15s.
% g_burgersfun.m         - ADiMat generated forward mode file (with loops)
% g_burgersfun_noloop.m  - ADiMat generated forward mode file (without loops)