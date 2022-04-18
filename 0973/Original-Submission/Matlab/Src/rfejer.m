function [ x , Err , Rel, y] = rfejer( sgl, fun , varargin)

%RFEJER Semi-Automatic Rational Fejer Quadrature
%   X = RFEJER(SGL) returns the nodes X(1,:) and weights X(2,:) in the 
%   (N+M+1)-point rational interpolatory Fejer quadrature rule on the 
%   interval [-1,1] based on the sequence of poles SGL, where 
%   N = lenght(SGL) and M denotes the number of complex poles in SGL for 
%   which the complex conjugate was not included in SGL. (Hence, for each 
%   complex pole in SGL its complex conjugate will be inserted whenever the 
%   latter was not included in SGL.) The poles must be outside the interval 
%   [-1,1], may be infinite (which corresponds to polynomial interpolation 
%   upto degree DEGR == sum(SGL==inf) ), and complex poles do not need 
%   to be given in adjacent positions, as the order will be rearranged if
%   necessary. X(3,1:end-1) contains the actual sequence of poles that is 
%   used (i.e., SGL with inclusion of the possible missing complex  
%   conjugates, and with the performed rearrangements). 
%
%   X = RFEJER(SGL,FUN) tries to approximate the integral of a 
%   (matrix-valued) function FUN over the interval [-1,1] as accurate as 
%   possible within a limited number of iterations - bounded by the length
%   of SGL (after inclusion of possible missing complex conjugates) - over 
%   the number of nodes and weights in the rational interpolatory Fejer 
%   quadrature rule with poles among SGL. The function y=FUN(x) should 
%   accept a vector argument x. 
%
%   [X,ERR,REL] = RFEJER(...) also returns an accuracy estimate ERR and a 
%   parameter REL to indicate for each computed value in X why the 
%   iterations have stopped:
%       REL == 1 ==> The iterations stopped due to numerical convergence
%       REL == 2 ==> The iterations stopped due to the fact that the 
%                    maximal number of iterations have been reached
%       REL == 3 ==> The iterations stopped due to deterioration in the
%                    construction of the rational Fejer quadrature formulae 
%
%   [X,ERR,REL,Y] = RFEJER(...) also returns the nodes Y(1,:) and the
%   weights Y(2,:), together with the actual sequence of poles
%   Y(3,1:end-1), from the last iteration. In the case of a vector- or
%   matrix-valued function FUN, recomputing X by means of Y and
%   FUN may lead, however, to a result that differs from X due to detection
%   of numerical convergence (REL == 1) for one or more functions at an
%   earlier iteration.   
%
%   X = RFEJER(SGL,FUN,PARAM1,VAL1,PARAM2,VAL2,...) performs the
%   integration with specified values of optional parameters. The available
%   parameters are
% 
%   'Array', which takes the value TRUE or FALSE and should be used
%       whenever the function y=FUN(x) does not accept a vector argument x. 
%       By default ARRAY == FALSE. 
%   'Tol', relative error tolerance. RFEJER attempts to satisfy 
%       ERR <= K*Tol, where K denotes the number of nodes and weights in 
%       the rational interpolatory Fejer quadrature formula. By default 
%       TOL == 5 * eps.        
%   'Nmax', maximal number of nodes and weights that should be used in the
%       rational interpolatory Fejer quadrature formula to approximate the 
%       integral. By default NMAX == N+M+1. Whenever it is used, the 
%       maximal number will be set to NMAX == min{'Nmax'+c,N+M+1}, where c 
%       is either 0 or 1 (depending on whether SGL contains complex poles). 
%   'Nmin', minimal number of nodes and weights that should be used in the
%       rational interpolatory Fejer quadrature formula to approximate the 
%       integral. By default NMIN == 1. Whenever it is used, the minimal 
%       number will be set to NMIN == min{'Nmin','Nmax'-2-d,N+M-1-d}, where
%       d is either 0, 1 or 2 (depending on whether SGL contains complex 
%       poles) in order to ensure at least three iterations.
%
%   Remark:
%       Whenever FUN contains one or more constant functions, the 
%       parameter 'Array' should take the value TRUE.
%
%   Example 1: A scalar-valued function on the interval [-1,1]
%     omeg = 0.1;
%     sgl = sqrt(-1)*omeg*[1:1:16];
%     NumInt = rfejer(sgl,@(x)fx(x,omeg));
%   where the file fx.m defines the function:
%       %------------------------------%
%       function y = fx(x,w)
%       y = (pi*x/w)./sinh(pi*x/w);
%       %------------------------------%
%     
%   Example 2: A matrix-valued function on the interval [0,1]
%     fx = @(x) [1./(x+0.2),1./(x+0.2).^2; 1./(x+0.1).^3,1./(x.^2+0.01)];
%     sgl = -0.1*[1,1,1,sqrt(-1),2,2];
%     ival = [0,1];
%     [fx2,sgl2] = transf(fx,sgl,ival);
%     [NumInt,err] = rfejer(sgl2,fx2);
%
%   Example 3: Combining rational and polynomial interpolation
%     sgl = [0.2*sqrt(-1),inf*ones(1,10)];
%     [X,err] = rfejer(sgl);
%     fx = @(x) exp(x)./(25*x.^2+1);
%     NumInt = fx(X(1,:))*X(2,:)';
%     NumInt2 = rfejer(sgl,fx);
%
%   See also RCHEB, ERRW, TRANSF.

% -------------------------------------------------------------------------
% Authors: Karl Deckers & Ahlem Mougaida & Hedi Belhadjsalah
%
% Reference: Extended Rational Fejer Quadrature Rules based on 
%                                   Chebyshev Orthogonal Rational Functions
%
% Software revision date: August 25, 2016
% -------------------------------------------------------------------------

% constants
pinf = 5 / eps;	  % poles >= pinf considered equal to inf
ptol = 5 * eps;   % poles closer than ptol apart considered equal
itolw = 100 * eps;  % tolerance on accuracy of the computed weights
itol = 5 * eps;   % tolerance on the computed approximations
ni = 1;           % initial iteration
ne = inf;         % last iteration

% check input
if nargin < 1,
    error('Not enough input arguments.');
else
    if any(isnan([real(sgl),imag(sgl)])),
        error('NaN is not a valid pole.');
    elseif any((abs(imag(sgl)) <= ptol) & (abs(real(sgl)) <= 1 + ptol)), 
        error('SGL contains poles too close to the interval IVAL.');
    else
        sgl = sgl(:).';
    end
    if nargin > 1,
        if ischar(fun),
            F = str2func(fun);
        else
            F = fcnchk(fun);
        end
        sizeF = size(F(0));
 
        % Process optional input.
        p = inputParser;
        p.addParamValue('Array',false,@validateArray);
        p.addParamValue('Tol',5*eps,@validateTol);
        p.addParamValue('Nmin',1,@validateNmin);
        p.addParamValue('Nmax',inf,@validateNmax);
        p.parse(varargin{:});
        optionStruct = p.Results;
        array = optionStruct.Array;
        itol = optionStruct.Tol;
        ni = optionStruct.Nmin;
        ne = optionStruct.Nmax;
        
        if ~array,
            F([-0.1,0.1]);
        end
    end
end

% process poles

sgl=pcheck(sgl,ptol);
sgl(abs(sgl) >= pinf) = inf;
mint = length(sgl);
b = real(sgl)+sqrt(-1)*abs(imag(sgl));
b = 1 ./ (b + (2 * (real(b)>=0) - 1) .* sqrt(b.^2 - 1));
index = find(imag(sgl)<0);
b(index) = conj(b(index));

% minimal and maximal number of points

if (mint < 3) || ((mint == 3) && ~(isreal(sgl(end-1)))),
    ni = 1; ne = mint;
else
    index = 1+sum((imag(sgl(1:3)) ~= 0));
    ne = min([max([ne-1,index,2]),mint]);
    if (ne < 3),
        ni = 1;
    else
        index = (imag(sgl(ne-2:ne)) ~= 0);
        index = ne - max([1,sum(index)]); % we need at least 3 iterations
        ni = min([ni,index]); 
    end
end
   
% initialise memory

vk = zeros(mint+1,1); % Kinv * Re{J(\varphi)}
vk(1) = 2/sqrt(pi);
k = 1;
if nargin > 1,
    x = zeros(sizeF);
    FQ = 2*F(0); % 1-point quadrature
else
    x = zeros(3,mint+1);
end
Err = inf*ones(size(x)); 
Rel = zeros(size(x)); 
conv = min(Rel(:));
war = 0;
p1 = zeros(size(x));
p2 = p1;
mint2 = ne;

% compute rational Chebyshev nodes and weights

% warning('off', 'Octave:possible-matlab-short-circuit-operator');
[xM,LM,erc]=rcheb([sgl,inf]);
% warning('on', 'Octave:possible-matlab-short-circuit-operator');
RphiM = mQ(b,xM); % compute Re{Q}^T

% iterations
while  conv==false &&  k<(mint2+1), 
    ak = sgl(k);
    rc = 2-isreal(ak);   % real/complex
    
    % initial step for new pole
    mk=length(find(ak==sgl(1:k)));
    if mk==1, 
        pos=find(ak==sgl(k:end)); % multiplicity
        m=length(pos);
        if rc==2,
            posc = pos;
            pos = zeros(1,2*m);
            pos(1:2:end) = posc;
            pos(2:2:end) = find(conj(ak)==sgl(k:end)); 
        end
        vk(pos+k) = JnB(ak,m,rc,ptol); % integral non-orthogonal basis
    end
 
    % compute Kinv * Re{J(\varphi)}

    vk(k+1:k+rc)= JoB(ak,mk,k+rc,xM,LM,vk(1:k+rc),RphiM,rc);
    k = k+rc;

    % approximate the integral
 
    if (nargin > 1) && ~(k < ni),
        if k < mint+1,
            % warning('off', 'Octave:possible-matlab-short-circuit-operator');
            [xi,lambda]=rcheb([sgl(1:k-1),inf]);
            % warning('on', 'Octave:possible-matlab-short-circuit-operator');
            Phi = mQ(b(1:k-1),xi);
        else
            xi = xM;
            lambda = LM;
            Phi = RphiM;
        end       
        W = lambda.*(Phi*vk(1:k))';
        war = 2*(abs(sum(W)-2)>k*itolw);
        if (war > 0) && (k-rc > ni),
            conv = true;
        else
            In1 = FQ;
            
            % (k+rc)-point quadrature
            FQ = evalquad(F,sizeF,xi,W,array);
            index = find(Rel==1);
            FQ(index)=In1(index);
            
            % convergence criterion
            In1 = abs(FQ-In1);
            aFQ = abs(FQ);
            index = find(aFQ > 0);
            In1(index) = In1(index)./aFQ(index);  
            index = find(Rel==1);
            In1(index) = Err(index);
            Err = In1;             
                % Put `%' at the beginning of line 259 and 260 to force the
                % iterations to continue until divergence has been detected
                % at line 238
            [Rel,p1,p2,mint2] = cvg(Err,p1,p2,k,itol,ne,mint2,Rel);
            conv = min(Rel(:));

            y = [xi;W;[sgl(1:k-1),inf]];
        end
    end
end

% output

if nargin == 1,
    x = [xM; LM.*(RphiM*vk)';[sgl,inf]];
    Err(1,:) = abs(erc);
    Rel(1,:)=1;
    Err(2,:)=abs(sum(x(2,:))-2);
    war = 2*(Err(2,1)>k*itolw);
    if ( ( war>0 ) || ( nargout>1 ) ),
        Err(2,:)=errW(x(3,:),x(1:2,:));
    end
    Rel(2,:)=1+war;
    Err(3,:)=0;
    y = x;
else
    x(1:end) = double(FQ);
    Err(1:end) = double(Err);
    index = find(Rel==0);
    if ~isempty(index),
        war = war + (k>mint2);
        Rel(index)=1+min([war,2]);
    end
    Rel(1:end) = double(Rel);
end

if nargout < 2,
    E = max(Err(:));
    if (war == 1) || ((war > 1) && (nargin == 1)),  
        warning('The desired accuracy of 1.e-13 may not be achieved.')
        warning('Estimated relative error on the approximation: %d',E)
    elseif (war > 1),        
        warning('The iterations failed to converge.')
        warning('Estimated relative error on the approximation: %d',E)
    end
end  

end

% -------------------------------------------------------------------------
function b=pcheck(a,tol)
%Check on complex conjugate pairs and on multiple poles

j=1;
m=length(a); 
while(j<=m)
    d = find(abs(1-a(j+1:m)/a(j))<tol);
    if ~isempty(d) 
            a(d(1)+j)=a(j);
    end
    j=j+1;
    if ~isreal(a(j-1))
       d=find(abs(1-a(j:m)/conj(a(j-1)))<tol);
       if ~isempty(d)
                a(d(1)+j-1)=[];
                m=m-1;
       end
       a(j:m+1)=[conj(a(j-1)),a(j:m)];
       m=m+1;
       j=j+1;
    end
end
b=a;

end

% -------------------------------------------------------------------------
function J = JnB(a,M,N,tol)
% Compute int_{-1}^1 [h_{a}(x)]^j dx for j=1,...,M 
%  Construct a column vector of length N*M where either 
%  N=1 (a is real) or N=2 (a is complex)  

c = abs(a)-1;
if isinf(a),
    J=zeros(1,M);
    J(2:2:end)=2./[3:2:M+1];

elseif (c<1e10*tol) || (log(10)*(4*M+50)+90*log(c)<0) || ((M<3) && (c<9)),
    J = (a^2-1)*log((a+1)/(a-1))-2*a; 
    if M>1, %forward recurrence
        F = 2;       
        for k=2:M,
            F1 = (a^2-1)*(1-(-1)^(k-1))/(k-1)-2*a*J(end)-a^2*F;
            F = J(end);
            J = [J;F1];
        end
    end 
else %series expansions
    f = mod(M,2);
    J=expansion(a,M,f,tol);
    if M>1, %backward recurrence
        F=J;
        J=expansion(a,M-1,1-f,tol);
        J=[J;F];
        for k=M-2:-1:1,
            F1=-((1/a)^2-1)*(1-(-1)^(k+1))/(k+1)-2*(1/a)*J(1)-(1/a)^2*F;
            F=J(1);
            J=[F1;J];
        end       
    end
end
 % complex case
if N == 2
     F = zeros(2*M,1);
     F(1:2:end) = J;
     F(2:2:end) = sqrt(-1)*J;   
     J = real(F);
end

end

%--------------------------------------------------------------------------
function J=expansion(a,M,f,tol)

j = 1;
while(log(abs(a)^(2*j)*(2*j+M+3-f)*(2*j+M+1-f)*tol)<0),
    j=j+1;
end
S = ([0:2:2*j]+M+3-f).*([0:2:2*j]+M+1-f).*a.^([0:2:2*j]);
S=sum(S.^(-1));
J = 2*(1-f)/(M+1)+(-1)^f*4*M/a^(2-f)*S; 

end

% -------------------------------------------------------------------------
function Q = mQ(beta,x)
% Construct matrix Re{Q}^T

n=length(x);
Q = zeros(n);
Q(1,:) = 1/sqrt(pi)* ones(1,n);
z = x+sqrt(-1)*sqrt(1-x.^2);
N=sqrt(1/(2*pi))*sqrt(1-abs(beta).^2);
B(1,:)=ones(1,n);
Bc(1,:)=B(1,:);
Q(2,:)=N(1)*(z.*Bc(1,:)./(1-beta(1).*z)+ 1./((z-beta(1)).*B(1,:)));
for k=2:n-1,
    B(k,:)=B(k-1,:).*(z-beta(k-1))./(1-conj(beta(k-1)).*z);
    Bc(k,:)=Bc(k-1,:).*(z-conj(beta(k-1)))./(1-beta(k-1).*z);
    Q(k+1,:)=N(k)*(z.*Bc(k,:)./(1-beta(k).*z)+ 1./((z-beta(k)).*B(k,:)));
end
Q=real(Q)';

end

% -------------------------------------------------------------------------
function J= JoB(a,m,k,x,lambda,Jf,Q,rc)
% Compute Kinv*int_{-1}^1 Re{phi_a(x)} dx 

if isinf(a),
    f = x.^m;
else
    f=((1-a*x)./(x-a)).^m;
    if rc==2,
        f=[real(f);-imag(f)];
    end
end

% Construct B
L=repmat(lambda,rc,1);
b = L.*f;
b = b*Q;
b = real(b(:,1:k));

%Compute J

beI=b(:,end-rc+1:end);
warning('off')
J = beI\(Jf(end+1-rc:end)-b(:,1:end-rc)*Jf(1:end-rc));
warning('on')

end

% -------------------------------------------------------------------------
function [Cn,p1,p2,mint2] = cvg(En,p1,p2,n,tol,mint,mint2,Co)
% Check convergence of approximation

  p1 = p2;
  p2 = Co;
  Cn = p2;
  index = find(En< sqrt(tol));  
  p2(index) = 1;
  index = find((p1+p2)==2);
  Cn(index) = (En(index) < n*tol);
  pM = min(p2(:));
  if  pM == 1,
      mint2 = min(max( n+8 , 2*n),mint);
  else
      mint2 = mint;
  end

end

% -------------------------------------------------------------------------
function FQ = evalquad(F,sizeF,xi,W,Array)
% Evaluate quadrature formula

n = length(xi);
if ~Array,
    mL = min([sizeF,n]);
    if mL==sizeF(2),
        G = F(xi);
        for j = 1:mL,
            FQ(:,j) = G(:,(j-1)*n+1:j*n)*W';
        end
    elseif mL==sizeF(1),
        G = F(xi');
        for j = 1:mL,
            FQ(j,:) = W*G((j-1)*n+1:j*n,:);
        end
    else
        Array = true;
    end
end
if Array,
    FQ = zeros(sizeF);
    for j = 1:n, 
        FQ = FQ+F(xi(j))*W(j);
    end
end

end

%--------------------------------------------------------------------------
function p = validateArray(x)
if ~(islogical(x))
    error(message('MATLAB:rfejer:invalidArray'));
end
p = true;
end

%--------------------------------------------------------------------------
function p = validateTol(x)
if ~(isfloat(x) && isscalar(x) && isreal(x) && x >= 0)
    error(message('MATLAB:rfejer:invalidTol'));
end
p = true;
end

%--------------------------------------------------------------------------
function p = validateNmin(x)
if ~(isfloat(x) && isscalar(x) && isreal(x) && x >= 0)
    error(message('MATLAB:rfejer:invalidNmin'));
end
p = true;
end

%--------------------------------------------------------------------------
function p = validateNmax(x)
if ~(isfloat(x) && isscalar(x) && isreal(x) && x >= 0)
    error(message('MATLAB:rfejer:invalidNmax'));
end
p = true;
end
