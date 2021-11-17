% problem.string.m
fprintf('\n====================================================================================\n')
fprintf('\n               >>>   Example 4a - Talbot Suite DE   <<<\n')
fprintf('\n         Application of Talbot''s method to solve the PDE equation\n')
fprintf('\n\t\t\tu_tt (x,t) = u_xx (x,t),         0 < x < L,  t>0\n')
fprintf('\n\twith conditions:')
fprintf('\n\t\t\tu(x,0+)   = x*sin(3*x)/6')
fprintf('\n\t\t\tu_t(x,0+) = sin(3*x)/6 + x*cos(3*x)/2')
fprintf('\n\t\t\tu(0,t)    = t*sin(3*t)/6\n')
fprintf('\n\tThe analytical solution is\n\n')
fprintf('\t\tu(x,t) = (x+t)*sin(3*(x+t))/6\n')
fprintf('\n\tand its Laplace Transform is\n\n')
fprintf('\t\t         sin(3*x)+3*x*cos(3*x)+s*x*sin(3*x)   3*sin(3*x)-s*cos(3*x)\n')
fprintf('\t\tU(x,s) = ---------------------------------- - ---------------------\n')
fprintf('\t\t                      6*(s^2+9)                     (s^2+9)^2\n')
fprintf('\n\twith double poles at +/-3i\n')
fprintf('\n\t     and abscissa of convergence: sigma0 = 0\n')
fprintf('\n====================================================================================\n')
