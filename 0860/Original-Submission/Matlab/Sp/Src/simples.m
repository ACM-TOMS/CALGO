function T = simples(k,d)
% T = simples(k,d) generates the k^d color schemes to
% be used in the d_dimensional simplex subdivision
% T is a cell array where each cell contains a matrix with
% the colors of the color scheme.
for n=1:(k^d)
  x = dec2base(n-1,k,d);
  cor=0;
  for i = 1:k
    CS(i,1)=cor;
    for j=1:d
      if x(d+1-j)==(i-1)
        cor=cor+1;
      end
      CS(i,j+1)=cor;
    end
  end
  T{n} = CS;
end
%------------------------------------------
function y = dec2base(x,b,n)
% y = dec2base(x,b,n) converts the decimal number x to
% y in base b with n digits:
% y(1)+y(2)*b+y(3)*b^2+...+y(n)*b^(n-1) = x
for i = 1:n
  d=fix(x/b);
  y(i)=x-d*b;
  x=d;
end
