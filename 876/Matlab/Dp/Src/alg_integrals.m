function ai = alg_integrals(c,alpha)
% For k = 0,1,2, evaluate the integrals
%                     1
%       I_k = Integral  u^k / |u-c|^alpha du
%                    -1 
% The exponent alpha satisfies 0 < alpha < 1. 
%
% This function is called by both Fie and ntrpFie.  It would
% be convenient to place it in a private directory.

if abs(c) > 10 % Taylor expansion for large values of |c|.
    an = zeros(1,14);
    an(1) = 1;
    for j = 2:14
        k = j-1;
        an(j) = an(j-1)*(-alpha-k+1)/k;
    end
    denom = 1:2:15;
    cf0 = an(1:2:13) ./ denom(1:7);
    cf1 = an(2:2:14) ./ denom(2:8);
    cf2 = an(1:2:13) ./ denom(2:8);

    u = 1/c^2;
    I_0 = cf0(7); I_1 = cf1(7); I_2 = cf2(7);
    for j = 6:-1:1
        I_0 = u*I_0 + cf0(j);
        I_1 = u*I_1 + cf1(j);
        I_2 = u*I_2 + cf2(j);
    end
    rca = 2/abs(c)^alpha;
    I_0 = rca*I_0;
    I_1 = -rca*I_1/c;
    I_2 = rca*I_2;
else
    if c < -1
        I_0 = ((1-c)^(1-alpha) - (-1-c)^(1-alpha))/(1-alpha);
        I_1 = ((1-c)^(2-alpha) - (-1-c)^(2-alpha))/(2-alpha) + c*I_0;
        I_2 = ((1-c)^(3-alpha) - (-1-c)^(3-alpha))/(3-alpha) + 2*c*I_1 - I_0*c^2;
    elseif c > 1
        I_0 = ((c+1)^(1-alpha) - (c-1)^(1-alpha))/(1-alpha);
        I_1 = ((c-1)^(2-alpha) - (c+1)^(2-alpha))/(2-alpha) + c*I_0;
        I_2 = ((c+1)^(3-alpha) - (c-1)^(3-alpha))/(3-alpha) + 2*c*I_1 - I_0*c^2;
    else  % -1 <= c <= 1
        I_0 = ((1-c)^(1-alpha) + (1+c)^(1-alpha))/(1-alpha);
        I_1 = ((1-c)^(2-alpha) - (1+c)^(2-alpha))/(2-alpha) + c*I_0;
        I_2 = ((1-c)^(3-alpha) + (1+c)^(3-alpha))/(3-alpha) + 2*c*I_1 - I_0*c^2;
    end
end
ai = [I_0,I_1,I_2];

end % alg_integrals
