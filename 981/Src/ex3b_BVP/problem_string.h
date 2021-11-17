char *str = "\n====================================================================================\n"
            "\n               >>>   Example 3b - Talbot Suite DE   <<<\n"
            "\n         Application of Talbot's method to solve the PDE equation\n"
            "\n\t\t\tu_t (x,t) = u_xx (x,t),         0 < x < L,  t>0\n"
            "\n\twith conditions:"
            "\n\t\t\tu(x,0+)   = x*(x-1)"
            "\n\t\t\tu(0,t) = 2*t"
            "\n\t\t\tu(L,t) = 2*t + L*(L-1)\n"
            "\n\tThe analytical solution is\n\n"
            "\t\tu(x,t) = 2*t + x*(x-1)\n"
            "\n\tand its Laplace Transform is\n\n"
            "\t\t          2    x*(x-1)\n"
            "\t\tU(x,s) = --- + -------\n"
            "\t\t         s^2      s\n"
            "\n\twith a double pole at 0\n"
            "\n\t     and abscissa of convergence: sigma0 = 0\n"
            "\n====================================================================================\n";
