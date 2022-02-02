function [NsingsTOT, Nsings, SINGS, MULT, sigma0] = LTsings2()
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                EXAMPLE 5

        Application of Talbot's method to solve
        the following PDE problem

                u_t  = u_xx  + u_yy,            0 < x,y < 1,   t>0
                u(x,y,0+) = x*(x-1) + y*(y-1)
                u(0,y,t) = u(1,y,t) = 4*t + y*(y-1)
                u(x,0,t) = u(x,1,t) = 4*t + x*(x-1)

        The analytical solution is

                u(x,y,t) = 4*t + x*(x-1) + y*(y-1)

        and Laplace Transform is

                          4    x*(x-1) + y*(y-1)
                U(x,s) = --- + -----------------
                         s^2           s

        where s = 0 is a double pole.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
    NsingsTOT = 1; % total number of singularities s_j
    Nsings = 1;    % number of singularities s_j with Im(s_j) ge 0
    sigma0 = 0;    % abscissa of convergence
    SINGS = 0;     % singularity
    MULT  = 2;     % polar multiplicity
end
