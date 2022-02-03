%
%   File:              cg.m
%
%   Description:       Conjugate Gradient Method
%                      to solve a positive definite system
%                      of equations  Ax = b
%
%   Authors:           Marielba Rojas
%                        mr@imm.dtu.dk
%                      Sandra A. Santos
%                      Danny C. Sorensen
%
%   Version:           1.2
%
%   System:            MATLAB 6.0 or higher 
%
%   Date:              15 March 2007
%
%   Functions called:  matvec, norm, length, zeros
%
%   Call: x = cg(A,b,epsilon);
%

function x = cg(A,b,epsilon)

    n   = length(b);

    x   = zeros(n,1);
    rc  = b;
    p   = rc;
%
%   Main loop
%
    while ( norm(rc) > epsilon )

       v     = matvec(p,A);
       alpha = (rc'*rc)/(p'*v);
       x     = x + alpha*p;
       ra    = rc;
       rc    = rc - alpha*v;
       beta  = (rc'*rc)/(ra'*ra);
       p     = rc + beta*p;

    end
