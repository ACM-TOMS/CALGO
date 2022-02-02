% MONOMIAL2BERNSTEIN    monomial2Bernstein(c) computes the coefficients of 
% the polynomial p(t)=c(1)+c(2)t+...+c(n)t^(n-1) respect to the Bernstein 
% basis of n-1 degree
function bernsteinCoef = monomial2Bernstein(varargin)
if nargin==0 || nargin>2
    disp('??? Error using ==> monomial2Bernstein');
    disp('Wrong dimesions of the argument');
    return;
else
    monomialCoef = varargin{1};
    [m n] = size(monomialCoef);
end
switch nargin
    case 1
        if m==1 && n~=1
            degree = n-1;
            dim = 2;
        elseif m~=1 && n==1
            degree = m-1;
            dim = 1;
        elseif m==1 && n==1
            bernsteinCoef = monomialCoef;
            return;
        else
            disp('??? Error using ==> monomial2Bernstein');
            disp('Wrong dimesions of the argument');
            return;
        end
    case 2
        dim = varargin{2};
        if dim == 1 % each column of coefMonomial are the coefficients of a polynomial in power form
            degree = m-1;
        elseif dim == 2 % each row of coefMonomial are the coefficients of a polynomial in power form
            degree = n-1;
        else
            disp('??? Error using ==> monomial2Bernstein');
            disp('Wrong dimesions of the argument');
            return;
        end
end

bernsteinCoef = zeros(m,n);
switch dim
    case 1
        for i=0:degree
            bernsteinCoef(i+1,:) = monomialCoef(i+1,:)/nchoosek(degree,i);
        end
        for r=1:degree
            for i=degree+1:-1:r+1
                bernsteinCoef(i,:) = bernsteinCoef(i-1,:)+bernsteinCoef(i,:);
            end
        end
    case 2
        for j=0:degree
            bernsteinCoef(:,j+1) = monomialCoef(:,j+1)/nchoosek(degree,j);
        end
        for r=1:degree
            for j=degree+1:-1:r+1
                bernsteinCoef(:,j) = bernsteinCoef(:,j-1)+bernsteinCoef(:,j);
            end
        end
end
end

