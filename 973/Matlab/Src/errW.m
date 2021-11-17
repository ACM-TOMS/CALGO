function [ err ] = errW( sgl , x )

%ERRW Estimate for the error on the computed weights x(2,:) 
%   * SGL is a sequence of poles plus an additional pole at infinity
%   * x(1,:) and x(2,:) are the nodes and weights in the corresponding
%       quadrature formula.
%   * For each complex pole in SGL, its complex conjugate should be 
%       included too in SGL; hence, it should hold that
%       N == length(SGL) == lenght(x(i,:))  i = 1, 2 
%   * ERR is a vector of length N containing the estimated (relative) errors

n = length(sgl);
X = x(1,:);
L = x(2,:)';
F = zeros(n);
Jf = zeros(n,1);
F(1,:) = 1;
Jf(1)=2;
for j=1:n
    a = sgl(j);
    m = length(find(a==sgl(1:j)));
    if isinf(a)
        F(j,:)=x(1,:).^(m-1);
        Jf(j)=(1-(-1)^(m))/(m);
    else
        F(j,:)=1./(x(1,:)-a).^m;
        if m==1,
            Jf(j)=log(1-a)-log(-1-a);
        else
            Jf(j)=( (1-a)^(1-m) - (-1-a)^(1-m) )/(1-m);
        end
    end
end

err = abs(Jf - F*L);
index = find(~(Jf==0));
err(index)=err(index)./abs(Jf(index));
err = err';

end
