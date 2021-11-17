% problem.string.m
fprintf('\n====================================================================================\n')
fprintf('\n               >>>   Example 3a - Talbot Suite DE   <<<\n')
fprintf('\n         Application of Talbot''s method to solve the PDE equation\n')
fprintf('\n\t\t\tu_t (x,t) = u_xx (x,t),         0 < x < L,  t>0\n')
fprintf('\n\twith conditions:')
fprintf('\n\t\t\tu(x,0+)  = x*(x-1)')
fprintf('\n\t\t\tu(0,t)   = 2*t')
fprintf('\n\t\t\tu_x(0,t) = -1\n')
fprintf('\n\tThe analytical solution is\n\n')
fprintf('\t\tu(x,t) = 2*t + x*(x-1)\n')
fprintf('\n\tand its Laplace Transform is\n\n')
fprintf('\t\t          2    x*(x-1)\n')
fprintf('\t\tU(x,s) = --- + -------\n')
fprintf('\t\t         s^2      s\n')
fprintf('\n\twith a double pole at 0\n')
fprintf('\n\t     and abscissa of convergence: sigma0 = 0\n')
fprintf('\n====================================================================================\n')
