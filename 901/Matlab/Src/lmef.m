function [c,ct,Res]=lmef(e_i,a,n,m)
% Function lmef calculates the coefficients of
% a linear multistep method with exponential fitting.
% Usage:
% [c,ct,e]=lmef(e_i,a,[ n ],[ m ])
% INPUT:
%   e_i     : 'e' for explicit methods, 'i' for implicit methods and
%             'b' for backword differentiation formulas (BDF).  
%   a       : a column vector with the values of the 
%             a-coefficients if explicit or implicit method is
%             selected. In the case of BDF, the values of the 
%             vector are neglected and only the dimension is
%             taken into account.
%   n       : a row vector with length equal to the number of
%             frequencies to be fitted. A value 1, means the
%             exponential exp(qx) is fitted and a value of 2, means
%             both exp(qx) and exp(-qx) are fitted.
%   m       : a row vector with the levels of tuning of the given
%             frequencies.
% OUTPUT:
%   c       : The a- and b-coefficients of the method in a compact 
%             form c=(a;b);
%   ct      : The Taylor expansion of c-coefficients.
%   e       : The reminder of the method.
%
%   NOTE:   
%   1.  If n,m are missing, algebraic fitting is assumed. If m
%       is missing, all levels of tuning assumed to be 0.
%   2.  In case of exponential fitting, the values of c and ct
%       are symbolic expressions, in which 'h' stands for the
%       step size and 'q_1','q_2',...,'q_n' are the n-different
%       frequencies that are fitted.
%   3.  The residual of the method is a symbolic expression in
%       which 'h' stands for the step size of the method, 
%       'q_1','q_2',...,'q_n' are the n-different frequencies 
%       that are fitted and 'y_1','y_2',...,'y_n' are the 
%       derivatives of the estimated function 'y'.
%
%   USAGE:
%   1.  Explicit method:
%       [c,ct,e]=lmef('e',[1,-1,1,-1]');
%   2.  Implicit method:
%       [c,ct,e]=lmef('i',[1,-1,1,-1]');
%   3.  Explicit with trigonometric fitting:
%       [c,ct,e]=lmef('e',[1,-1,0,1,-1]',[2]);
%   4.  BDF with 5 steps:
%       [c,ct,e]=lmef('b',[1,1,1,1,1]');

%--------------------------------------------------------------------
% Customization parameters
%--------------------------------------------------------------------
    te_order=10;
    rc_order=1E-5;
    global h;
    h=sym('h','real');

%--------------------------------------------------------------------
% Check input parameters
%--------------------------------------------------------------------
if ((nargin < 2) | (nargin > 4))   % Wrong
    disp('********************************************************');
    disp('LMEF: Linear Multistep Methods with Exponential Fitting');
    disp('Usage:');
    disp('      [c,ct,e]=lmef(e_i , a , n , m);');
    disp('IN:');
    disp('  e_i :   ''e''-explicit  / ''i''-implicit / ''b''-BDF' );
    disp('  a   :   the coefficients - a');
    disp('  n   :   frequency vector');
    disp('  m   :   level of tuning vector');
    disp('********************************************************');
    error 'Wrong argument number.'
end
if nargin == 2  % Algebraic fitting only
    n=[1];
    m=[0];
end
if nargin == 3  % Exp fitting (all levels of tuning are zero)
    m=ones(1,max(size(n)));
end
if nargin == 4
    for j=1:max(size(m))
        m(j)=m(j)+1;
    end
end
if ((e_i~='e')&(e_i~='i')&(e_i~='b'))
    error('Wrong 1st input value. Must be ''e'' , ''i'' or ''b''.');
end
if max(size(n))~=max(size(m))
    error('The dimension of n and m must be the same');
end
for j=1:max(size(n))
    if ((n(j)~=1) & (n(j)~=2))
        error('Freq vector must contain 1''s and 2''s.');
    end
end
%--------------------------------------------------------------------
% Check the a-coefficients
%--------------------------------------------------------------------
k=max(size(a))-1;%This is the dimension of the problem
if (e_i~='b')
    disp('Checking the a-coefficients...');
    if (sum(a)~=0)
        error('Characteristic polynomial must have 1 as root.');
    end
    ra=roots(a);
    for j=1:k
        if (abs(ra(j)) > 1+rc_order)
            error('Characteristic polynomial roots must lie on unit disc.');
        end
        for jj=j+1:k
            if ( (ra(j)==ra(jj)) & (abs(ra(j))==1))
                error('Roots lying on unit circle must have multiplicity 1.');
            end
        end
    end
    disp('Ok');
end
%--------------------------------------------------------------------
% Build the conditions and solve for unknowns
%--------------------------------------------------------------------
disp('Calculating the coefficients...');
[C,ma]=cr_con(e_i,k,n,m);
if (ma<0)
    error('Incompatible system.');
end
c=solve_c(C,a,e_i);
a=c(1:k+1);
b=c(k+2:2*k+2);
disp('Ok');
%--------------------------------------------------------------------
% In case of exp-fitting, calculate Taylor expansions
%--------------------------------------------------------------------
ct=c;
if (n*m' > 0)
    disp('Calculating the Taylor expansion of coefficients...');
    syms ctc;
    if (e_i=='i')
        sj=k+2;
        fj=2*k+2;
    elseif (e_i=='e')
        sj=k+3;
        fj=2*k+2;
    else
        sj=1;
        fj=k+1;
    end
    for j=sj:fj
        jj=10;
        deg=0;
        while ((deg < te_order) & (ct(j)~=0))
            [ct(j),e]=my_mtaylor(c(j),jj,h);
            if (ct(j)==0)
                deg=0;
            else
                deg=eval(maple('degree',ct(j),h));
            end
            jj=jj+(te_order-deg);
        end
    end
    disp('Ok');
end
%--------------------------------------------------------------------
% Calculate the expression for residual
%--------------------------------------------------------------------
disp('Calculating the Residual...');
if (mod(k,2)==1)
    kk=k+2;
else
    kk=k+3;
end
TC=zeros(1,kk);
TC=sym(TC);
yt=zeros(1,kk);
yt=sym(yt);
for j=1:kk
    buf=sprintf('y_%d',j);
    yt(j)=sym(buf);
end
for j=1:k
    TC(1)=TC(1)+(-j)*a(j+1);
end
TC(1)=TC(1)-sum(b);
TC(1)=TC(1)*yt(1)*h;
for nn=2:kk
    s1=0;
    for j=1:k
        s1=s1+a(j+1)*(-j)^nn;
    end
    s2=0;
    for j=1:k
        s2=s2+nn*b(j+1)*(-j)^(nn-1);
    end
    TC(nn)=(s1-s2)*yt(nn)*h^nn/factorial(nn);
end
disp('Ok');
%--------------------------------------------------------------------
% Calculate the Taylor expansion of residual
%--------------------------------------------------------------------
disp('Calculating the Taylor expansion of the Residual...');
TCC=sum(TC);
Res=0;
j=1;
e=1;
while ((Res==0) & (e==1))
    [Res,e]=my_mtaylor(TCC,j,h);
    j=j+1;
end
if e==1
    TCC=Res;
    Res=0;
    j=1;
    e=1;
    while ((Res==0) & (e==1))
        [Res,e]=my_mtaylor(TCC,j,h);
        j=j+1;
    end
    Res=taylor(TCC,j,h);
    Res=simplify(Res);
else
    disp('Cannot calculate Reminder.');
end
disp('Ok');
return
%********************************************************************
%********************************************************************
%********************************************************************
function [C,ma]=cr_con(e_i,k,n,m)
%-------------------------------------------------
%--Function cr_con
% This function calculates the conditions that 
% the coefs a's and b's must satisfy in a Linear
% Multistep method with exponential fitting.
% The condition À»a_j=0 is assumed to hold is assumed
% in the case of explicit or implicit method.
% The meaning of the matrix C is that C*(a;b)=0.
% 
%--Input:
%   e_i :   explicit, implicit or BDF method
%   k   :   order of the method. This means that
%           we have (k+1) equations to build.
%   n   :   number of different exponentials
%           to fit.
%   m   :   vector containing the multiplicities
%           of the n different exponents.
%--Output:
%   C   :   (k+1) x (2*k+2) array containing the
%           conditions for the coef's            
%   ma  :   the order of algebraic fitting
%-------------------------------------------------
global h;
% C is the output matrix
C=zeros(k+1,2*k+2);
C=sym(C);
% Check first for system compatibility
exp_eqs=n*m';   % This is the total number of
                % equations resulting from
                % exponential and trigonometric
                % fitting.
if ((e_i=='i'))
    ma=k+1-exp_eqs; % ma is the total number of
else                % algebraic conditions
    ma=k-exp_eqs;   % which is the number of 
end                 % unknowns minus exp_eqs.
if (ma<0)       % Incompatible system
    return;
end
%
% Definition of algebraic initial condition
Al=zeros(1,2*k+2);
for j=1:k+1
    Al(j)=1;
end
%
% Condition for explicit method
if e_i=='e'
    Ae=zeros(1,2*k+2);
    Ae(k+2)=1;
end
%
% Algebraic conditions
for j=1:ma
    C(j,:)=Al*Lambda(k,j);
end
off=ma+1;
%
% Explicit condition
if e_i=='e'
    C(off,:)=Ae;
    off=off+1;
end
%
% BDF Formula
if e_i=='b'
    C(off,:)=Al;
    off=off+1;
end
%
% Conditions for exponential fitting
w=sym('w');
q=sym('q');% q is the frequency
%
% Exponential condition
Aex=zeros(1,2*k+2);
Aex=sym(Aex);
for j=1:k+1
    Aex(j)=exp((k/2-j+1)*q*h);
    Aex(j+k+1)=-q*h*exp((k/2-j+1)*q*h);
end
%
% Trigonometric condition
Atr=zeros(2,2*k+2);
Atr=sym(Atr);
for j=1:k+1
    Atr(1,j)=Aex(j);%sin((k/2-j+1)*u);
    Atr(1,j+k+1)=Aex(j+k+1);%-u*cos((k/2-j+1)*u);
    Atr(2,j)=subs(Aex(j),q,-q);%cos((k/2-j+1)*u);
    Atr(2,j+k+1)=subs(Aex(j+k+1),q,-q);%u*sin((k/2-j+1)*u);
end
%
% Build exp/trig conditions
for j=1:max(size(n))
    if (n(j)==1)    % Exponential fitting
        buf=sprintf('q_%d',j);
        q_d=sym(buf,'unreal');
        Tex=subs(Aex,q,q_d);
        for jj=0:m(j)-1
            C(off,:)=Tex*Lambda(k,jj);
            off=off+1;
        end
    else            % Trigonometric fitting
        buf=sprintf('q_%d',j);
        q_d=sym(buf,'unreal');
        Ttr=subs(Atr,q,q_d);
        for jj=0:m(j)-1
            C(off:off+1,:)=Ttr*Lambda(k,jj);
            off=off+2;
        end
    end
end
return
%********************************************************************
%********************************************************************
%********************************************************************
function N=N_k(order,expon)
%-------------------------------------------------
%--Input:
%   order   :   order of the correspondin method
%   expon   :   the power of the natrix
%--Output:
%   N       :   (order+1) x (order+1) diagonal 
%               array with N(j,j)=(j-1)^expon
%-------------------------------------------------
    N=zeros(order+1,order+1);
    if expon<=0
        N=eye(order+1);
    else
        for j=1:order+1
            N(j,j)=(j-1)^expon;
        end
    end
return
%********************************************************************
%********************************************************************
%********************************************************************
function L=Lambda(order,expon)
%-------------------------------------------------
%--Input:
%   order   :   order of the correspondin method
%   expon   :   the power of the natrix
%--Output:
%   L       :   (2*order+2) x (2*order+2) 
%               triangular array equal to
%  [N_k^expon , expon*N_k^(expon-1);0 , N_k^expon]
%-------------------------------------------------
    if expon<=0
        L=eye(2*order+2,2*order+2);
    else
        L=zeros(2*order+2,2*order+2);
        L(1:order+1,1:order+1)=N_k(order,expon);
        L(1:order+1,order+2:2*order+2)=expon*N_k(order,expon-1);
        L(order+2:2*order+2,order+2:2*order+2)=N_k(order,expon);
    end
return
%********************************************************************
%********************************************************************
%********************************************************************
function c=solve_c(C,a,e_i)
%-------------------------------------------------
%--Input:
%   C       :   (k+1)x(2k+2) matrix such that
%               C*(a;b)=0
%   a       :   the coefs a's
%   e_i     :   the type of the method
%--Output:
%   c       :   the calculated coefficients
%-------------------------------------------------
%
% Check first if matrix C is OK
[k,kk]=size(C);
if 2*k~=kk
    error('Wrong matrix C size');
end
ka=max(size(a));
if ka~=k
    error('Wrong vector a size');
end
%
% Split C into (C1:C2);
C1=C(1:k,1:k);
C2=C(1:k,k+1:kk);
%
% If implicit or explicit, solve for b's
if (e_i~='b')
    F=-C1*a;
    b=C2\F;
else
%
% If BDF, solve for a's. b=[1 0 ...]
    b=zeros(k,1);
    b(1)=1;
    F=-C2*b;
    a=C1\F;
end
c=[a;b];
return
%********************************************************************
%********************************************************************
%********************************************************************
function [bt,e]=my_mtaylor(b,m,h)
%-------------------------------------------------
%--Input:
%   b       :   an expression of h
%   m       :   the order of the expansion
%   h       :   the free variable
%--Output:
%   bt      :   the Taylor expansion of b
%   e       :   1 if the expansion exists
%               0 otherwise
%-------------------------------------------------
syms gn gp gnt gpt;
e=1;bt=0;
%
%Split given expression into numerator and denominator
[gn,gp]=numden(b);
%
%Calculate Taylor expansion of numerator
l=limit(gn,h,0);
gnt=l;
for j=1:m
    aux=(gn-gnt)/h^j;
    l=limit(aux,h,0);
    gnt=gnt+l*h^j;
end
%
%Calculate Taylor expansion of denominator
l=limit(gp,h,0);
gpt=l;
for j=1:m
    aux=(gp-gpt)/h^j;
    l=limit(aux,h,0);
    gpt=gpt+l*h^j;
end
%
%Check if expansion does not exists
if gpt==0
    if gnt==0
        e=1;
    else
        e=0;
    end
    return;
end
%
%Finally calculate Taylor expansion of fraction
aux=gnt/gpt;
l=limit(aux,h,0);
bt=l;
for j=1:m
    aux1=(aux-bt)/h^j;
    l=limit(aux1,h,0);
    bt=bt+l*h^j;
end
return
