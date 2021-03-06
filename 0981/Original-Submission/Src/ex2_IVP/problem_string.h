char *str = "\n====================================================================================\n"
            "\n               >>>   Example 2 - Talbot Suite DE   <<<\n"
            "\n         Application of Talbot's method to solve the PDE equation\n"
            "\n\t\t\tu_tx (x,t) = exp(-x)*cos(t),         x>0,  t>0\n"
            "\n\twith conditions:"
            "\n\t\t\tu_x(x,0+) = 0"
            "\n\t\t\tu(0,t)    = 0\n"
            "\n\tThe analytical solution is\n\n"
            "\t\tu(x,t) = sin(t)*(1 - exp(-x))\n"
            "\n\tand its Laplace Transform is\n\n"
            "\t\t         1 - exp(-x)\n"
            "\t\tU(x,s) = -----------\n"
            "\t\t           s^2 + 1\n"
            "\n\twith simple poles at +/-i\n"
            "\n\t     and abscissa of convergence: sigma0 = 0\n"
            "\n====================================================================================\n";
