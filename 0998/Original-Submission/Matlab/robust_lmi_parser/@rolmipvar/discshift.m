function polyout = discshift(varargin)
% [polyout] = discshift(polyin,eta,bounds,targetin,targetout)   
% Applies up to eta time-shifts in a discrete-time parameter-varying
% polynomial polyin, contained in the simplex domaing targetin. The input
% bounds is a Nx1 vector that informs the bounds for the variation rates of
% all the N parameters, i.e., dm_i <= a_i(k+1) - a_i(k) <= dM_i. The resultant
% polynomial polyout is a cell structure with eta+1 positions, where
% polyout{n} corresponds to the input polynomial dependent on a(k+n-1) and
% is in the new simplex domain targetout.
%
% The syntax
% [polyout] = discshift(polyin,eta)
% may also be used if the input polynomial depends on only one simplex
% domain, with parameters varying under arbitrary variation rates.
%
% [polyout] = discshift(polyin,eta,bounds)
% This syntax may be used if the input polynomial depends on only one
% simplex domain, but the variation rates of the parameters are bounded.
% If, for instance, dm_1 <= a_1(k+1) - a_1(k) <= dM_1 and
% dm_2 <= a_2(k+1) - a_2(k) <= dM_2, then bounds = [dm_1 dM_1; dm_2 dM_2]
% 
%
% [polyout] = discshift(polyin,eta,bounds,targetin)
% The resultant simplex will be allocated to the first available simplex
% index.
% 
% [polyout] = discshift(polyin,eta,targetin)
% Same as the prior syntax but considering arbitrary variation rates.
%
% [polyout] = discshift(polyin,eta,targetin,targetout)
% This syntax may be used if the variation rates of the parameters are
% arbitrary.
%
% Important notice: The command discshift needs the MPT3 package to be
% installed (https://www.mpt3.org/)

if (~exist('Polyhedron'))
    error('The command discshift needs the MPT3 package to be installed (https://www.mpt3.org/)');
end

if ((nargin < 2) || (nargin > 5))
    error('Input error. Type ''help discshift'' for more information');
else
    polyin = varargin{1};
    eta = varargin{2};
    if (~isa(polyin,'rolmipvar'))
        error('The first input must be a rolmipvar variable');
    end
    
    if (nargin == 2)
        targetin = 1;
        bounds = [-ones(polyin.vertices(targetin),1) ones(polyin.vertices(targetin),1)];
        targetout = length(polyin.data(1).exponent)+1;
    end
    
    if (nargin == 3)
        if (length(varargin{3}) > 1)
            bounds = varargin{3};
            if (length(polyin.data(1).exponent) == 1)
                targetin = 1;
                targetout = 2;
            else
                error('If the polynomial depends on more than one simplex domain, please inform the simplex to perform the shift operation');
            end
        else
            targetin = varargin{3};
            bounds = [-ones(polyin.vertices(targetin),1) ones(polyin.vertices(targetin),1)];
            targetout = length(polyin.data(1).exponent)+1;
        end
    end
    
    if (nargin >= 4)
        if (length(varargin{3}) > 1)
            bounds = varargin{3};
            targetin = varargin{4};
            if (nargin == 4)
                targetout = length(polyin.data(1).exponent)+1;
            else
                targetout = varargin{5};
            end
        else
            targetin = varargin{3};
            bounds = [-ones(polyin.vertices(targetin),1) ones(polyin.vertices(targetin),1)];
            targetout = varargin{4};
        end 
    end
end

N = polyin.vertices(targetin);

Abar = [kron(eye(eta+1),[1;-1]); [kron(eye(eta),[-1;1]) zeros(2*eta,1)]+[zeros(2*eta,1) kron(eye(eta),[1;-1])]];
A = [];
b = [];
for cont=1:N
    A = blkdiag(A,Abar);
    bar{cont} = [kron(ones(eta+1,1),[1;0]); [kron(ones(eta,1),[bounds(cont,2);-bounds(cont,1)])]];
    b = [b; bar{cont}];
end
Ae = kron(ones(1,N),eye(eta+1));
be = ones(eta+1,1);


P = Polyhedron('A',A,'b',b,'Ae',Ae,'be',be);
Haux = P.V';  %[a1(k); a1(k+1); a2(k); a2(k+1)]

em = [1 zeros(1,N*(eta+1)-1)];
M = zeros(N*(eta+1),N*(eta+1));
cont = 1;
for conteta=1:eta+1
    em = circshift(em,[0 conteta-1]);
    for contN=1:N
        M(cont,:) = em;
        em = circshift(em,[0 eta+1]);
        cont = cont + 1;
    end
    em = [1 zeros(1,N*(eta+1)-1)];
end
H = M*Haux; %[a1(k); a2(k); a1(k+1); a2(k+1)]

clear M;
for conteta = 1:eta+1
    M = H((conteta-1)*N + 1:conteta*N,:);
    polyout{conteta} = simpltransf(polyin,targetin,targetout,M);
    if (conteta > 1)
        polyout{conteta}.label = strcat(strcat(strcat(polyout{conteta}.label,'(k+'),num2str(conteta-1)),')');
    end
end

return