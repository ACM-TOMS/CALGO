function li = log_integrals(c)
% For k = 0,1,2, evaluate the integrals
%                     1
%       I_k = Integral (u^k)log|u-c| du
%                    -1 
%
% This function is called by both Fie and ntrpFie.  It would
% be convenient to place it in a private directory.

ua = abs(1+c); ub = abs(1-c); 
absc = abs(c); 

if absc > 10  % Taylor expansion for large values of |c|.
    logabsc = log(absc);
    u = 1/c^2;
    I_0 = 2*logabsc - u*(1/3 + u*(0.1 + u*(1/21 +u*(1/36 + u*(1/55 + u/78)))));
    I_1 = (-2/c)*(1/3 + u*(1/15 + u*(1/35 + u*(1/63 + u*(1/99 + u/143)))));
    I_2 = logabsc/1.5 - u*(0.2 + u*(1/14 + u*(1/27 + u*(1/44 + u*(1/65 + u/90)))));
elseif c == +1
    I_0 = 2*log(2) - 2;
    I_1 = -1;
    I_2 = -8/9 + log(2)/1.5;
elseif c == -1
    I_0 = 2*log(2) - 2;
    I_1 = 1;
    I_2 = -8/9 + log(2)/1.5;
else
    % Use this for smaller values of |c|, except c = +1,-1. As |c|
    % increases, the values are increasingly inaccurate, especially 
    % the values for I_1 and I_2.
    logua = log(ua); logub = log(ub);
    I_0 = (1-c)*logub + (1+c)*logua -2;
    I_1 = -c + 0.5*(1-c^2)*log(ub/ua);
    I_2 = -(1+3*c^2)/4.5 + ((1-c^3)*logub + (1+c^3)*logua)/3;
end
li = [I_0,I_1,I_2];

end % log_integrals