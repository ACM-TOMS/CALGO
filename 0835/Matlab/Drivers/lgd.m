function p = lgd(n)
%
%  test polynomial suggested by Goedecker
%   Legendre polynomial of degree n
%
    p0 = [1];
    p1 = [1,0];
    if n == 0
        p = p0;
    elseif n == 1
        p = p1;;
    elseif n > 1
        for k = 1:n-1
            q1 = (2*k+1)*[p1,0];
            m = length(q1)-length(p0);
            q0 = k*[zeros(1,m),p0];
            p = (q1-q0)/(k+1);
            %
            p0 = p1;
            p1 = p;
        end;
    end;
