function [Q,nfun]=da0glob(f,a,tolabs,tolrel,trace)
%DA0GLOB  Numerically evaluates an integral using adaptively
%   a Newton-Cotes 5/9 point rule. The routine is globally adaptive and the
%   fact that the routine makes use of two rules implies that it is doubly adaptive.
%   Except for some minor changes this routine is identical to the routine
%   COTEGLOB published by this author in 2002.
%
%   [Q,NFUN]=DA0GLOB(F,A) approximates the sum of the integrals of F(X) from
%   A(i) to A(i+1), i=1,...,length(A)-1, to a relative precision 10^-3. Note that length(A)
%   has to be at least 2. Q is the estimate of the integral and NFUN is the number 
%   of function values used by the code. F is a function handle.
%   The function F must return a vector of output values if
%   given a vector of input values.
%
%   [Q,NFUN]=DA0GLOB(F,A,TOLABS) integrates to an absolute
%   error of TOLABS.
%
%   [Q,NFUN]=DA0GLOB(F,A,TOLABS,TOLREL) integrates to an absolute
%   error of max(TOLABS,EST*TOLREL), where EST is the current estimate of the
%   integral. The code also estimates the problems NOISE level
%   and it will not allow the effective absolute error to be below this
%   NOISE level: thus max([TOLABS,EST*TOLREL,NOISE])
%
%   [Q,NFUN]=DA0GLOB(F,A,TOLABS,TOLREL,TRACE) displays the left
%   and right end points of the current interval, the partial integral 
%   and the rule used (5 or 9).
%
%   The code also counts the number of function values and it may stop due
%   to the fact that the maximum number (10 000) has been reached.
%   Furthermore the code may stop processing an interval if the interval
%   size becomes too small giving a warning that this has happened. Finally
%   the code checks that the tolerance is not below the problem's noise
%   level. In such a case the stopping criterion is modified and a warning
%   given.
%
%   Modified by
%   Terje O. Espelid 24/09/06

%   References: Berntsen and Espelid, TOMS 1991
%               Espelid, in Espelid and Genz, 1992.
%               Espelid, Bit 2003
%               Espelid, Report 266 February 2004
%   Check the number of input variables:
  if(nargin<5), trace=0; end
  if(nargin<4), tolrel=0; end
  if(nargin<3), tolabs=0;tolrel=0.001;end
% The weights of the rules and nullrules
w1=[1.5555555555555556e-01   7.1111111111111114e-01   2.6666666666666666e-01   7.1111111111111114e-01   1.5555555555555556e-01];
nw1=[1.2710311885185779e-01  -5.0841247540743117e-01   7.6261871311114682e-01  -5.0841247540743117e-01   1.2710311885185779e-01
  -3.3628324334270127e-01   6.7256648668540253e-01   0.0000000000000000e+00  -6.7256648668540253e-01   3.3628324334270127e-01
   5.6842242780997809e-01  -2.8421121390498905e-01  -5.6842242780997809e-01  -2.8421121390498905e-01   5.6842242780997809e-01
  -6.7256648668540253e-01  -3.3628324334270127e-01   0.0000000000000000e+00   3.3628324334270127e-01   6.7256648668540253e-01];
w2=[6.9770723104056437e-02   4.1537918871252205e-01  -6.5467372134038804e-02   7.4045855379188708e-01  -3.2028218694885363e-01   7.4045855379188708e-01  -6.5467372134038804e-02   4.1537918871252205e-01   6.9770723104056437e-02];
nw2=[1.1018547692345270e-02  -8.8148381538762158e-02   3.0851933538566756e-01  -6.1703867077133512e-01   7.7129833846416884e-01  -6.1703867077133512e-01   3.0851933538566756e-01  -8.8148381538762158e-02   1.1018547692345270e-02
  -4.2674651711845396e-02   2.5604791027107238e-01  -5.9744512396583560e-01   5.9744512396583560e-01   0.0000000000000000e+00  -5.9744512396583560e-01   5.9744512396583560e-01  -2.5604791027107238e-01   4.2674651711845396e-02
   1.1236757938944257e-01  -4.7756221240513091e-01   6.1802168664193424e-01   2.8091894847360643e-02  -5.6183789694721287e-01   2.8091894847360643e-02   6.1802168664193424e-01  -4.7756221240513091e-01   1.1236757938944257e-01
  -2.3112700627433047e-01   6.3559926725440874e-01  -2.3112700627433047e-01  -5.2003576411724362e-01   0.0000000000000000e+00   5.2003576411724362e-01   2.3112700627433047e-01  -6.3559926725440874e-01   2.3112700627433047e-01
   3.9111964345081668e-01  -5.8667946517622505e-01  -3.0730829128278453e-01   2.5143405650409645e-01   5.0286811300819290e-01   2.5143405650409645e-01  -3.0730829128278453e-01  -5.8667946517622505e-01   3.9111964345081668e-01
  -5.5619114160254801e-01   2.7809557080127401e-01   5.1646320291665171e-01   3.5755144817306656e-01   0.0000000000000000e+00  -3.5755144817306656e-01  -5.1646320291665171e-01  -2.7809557080127401e-01   5.5619114160254801e-01
   6.6477556470172228e-01   1.6619389117543057e-01  -1.8993587562906350e-01  -4.0361373571175990e-01  -4.7483968907265878e-01  -4.0361373571175990e-01  -1.8993587562906350e-01   1.6619389117543057e-01   6.6477556470172228e-01
  -6.4550259924248821e-01  -4.8412694943186613e-01  -3.2275129962124410e-01  -1.6137564981062205e-01   0.0000000000000000e+00   1.6137564981062205e-01   3.2275129962124410e-01   4.8412694943186613e-01   6.4550259924248821e-01];
% w1 and nw1 are the weights and nullrules associated with the 5 point Newton-Cotes rule.
% w2 and nw2 are the weights and nullrules associated with the 9 point Newton-Cotes rule.
% The coordinates of the four new points in a five point approximation
  Z=[-3,-1,1,3]/4;
% Define the number of intervals given
  nint=length(a)-1;
  if nint < 1, warning('No interval is given.');Q=0;nfun=0;return,end
% Initialize: the function count, the maximum number of intervals, too small intervals
% (logical: interval) and the maximum number of function values allowed.
  nfun=0;int_max=2000;interval=0;nmax=10000;singular=0;
% Initialize several dynamic vectors for speed: The choice of value of
% int_max is connected to the fact that we have a bound 10 000 on the
% number of function values and if we use the 5 point rule everywhere
% (not very likely) then 2000 intervals times 4 f-evaluation per interval
% implies approx 8 000 function evaluations. Therefore, this will do in
% most situations.
  is=zeros(1,int_max);error=is;isa=is;lp=is;rp=is;flag=is;
  fu=zeros(int_max,9);
  for i=1:nint  
% Compute y and start the count of the number of function values
     x=linspace(a(i),a(i+1),9); y=feval(f,x);nfun=nfun+9;
% Check if some of the function values are infinite. In case there are: put
% the value(s) to zero.
     IN=find(isinf(y)); y(IN)=0;singular=max(1-isempty(IN),singular);
% Check if some of the function values are NaN: put
% those function values to 1 (any number might do).
     y(isnan(y))=1;
% Apply the degree 9 Cotes rule on the interval and
% initialize the datastructure by storing: error estimate, integral, absolute value, 
% function values, left and right point of the interval and
% signalling by a flag that a nine point rule has been used.
     [error(i),is(i),isa(i)]=cotes9(a(i),a(i+1),y,32,w2,nw2,interval);
     fu(i,1:9)=y;lp(i)=a(i);rp(i)=a(i+1);flag(i)=1;
  end
% Compute global values
  Q=sum(is);isabs=sum(isa);total=sum(error);noise=50*eps*isabs;

% global error test and check the number of function evaluations.

  while (total> max([abs(Q)*tolrel,tolabs,noise]))&&(nfun<nmax-3)

%  Find the interval with the greatest error estimate.
      [err,i]=max(error);
%  Prepare to process this interval: either by subdivision or by using a 9
%  point rule
      aa=lp(i);bb=rp(i);cc=(aa+bb)/2;total=total-err;Q=Q-is(i);isabs=isabs-isa(i);
      if flag(i)==1
%  On this interval a 9 point rule has been used. Subdivide it and use all nine function values
%  to compute: first the left part  of the integral
         z(1:5)=fu(i,1:5);[err,inte,inteabs]=cotes5(aa,cc,z(1:5),32,w1,nw1);zz(1:5)=fu(i,5:9);
         error(i)=err;is(i)=inte;isa(i)=inteabs;total=total+err;Q=Q+inte;isabs=isabs+inteabs;
         flag(i)=0;fu(i,1:2:9)=z(1:5);fu(i,2:2:8)=0;lp(i)=aa;rp(i)=cc;
%  and then the right part of the integral.
         nint=nint+1;i=nint;[err,inte,inteabs]=cotes5(cc,bb,zz(1:5),32,w1,nw1);
         error(i)=err;is(i)=inte;isa(i)=inteabs;total=total+err;Q=Q+inte;isabs=isabs+inteabs;
         flag(i)=0;fu(i,1:2:9)=zz(1:5);fu(i,2:2:8)=0;lp(i)=cc;rp(i)=bb;
      else
%  On this interval a 5 point rule has been used only. Compute f in 4 new
%  points and apply the 9 point rule
         hh=(bb-aa)/2;x4=cc+hh*Z;y4=feval(f,x4);nfun=nfun+4;
         IN=find(isinf(y)); y(IN)=0;singular=max(1-isempty(IN),singular);
         y(isnan(y))=1;
         fu(i,2:2:8)=y4(1:4);y(1:9)=fu(i,1:9);[err,inte,inteabs]=cotes9(aa,bb,y(1:9),32,w2,nw2,interval);
         error(i)=err;is(i)=inte;isa(i)=inteabs;flag(i)=1;
         total=total+err;Q=Q+inte;isabs=isabs+inteabs;
      end
% Redefine the global noise level
      noise=50*eps*isabs;
  end
  if (nfun>=nmax-3) 
       warning(['Stopping: maximum number of f-evaluations,'...
                'required tolerance may not be met.'])
  end
  if (interval==1)
         warning('Interval too small, required tolerance may not be met')
  end
  if max(abs(Q)*tolrel,tolabs) < noise
       warning(['Stopping: the tolerance is below the noise level of the problem,'...
                'required tolerance may not be met.']);
  end
  if singular==1
      warning('Singularity probably detected, required tolerance may not be met.')
  end
  if trace
     for i=1:nint
         disp([lp(i),rp(i),is(i),5+flag(i)*4])
     end
  end
  Q=sum(is);
  
function [err,Q,Qabs]=cotes9(a,b,y,D,w2,nw2,interval)
%
%   [err,Q,Qabs] = COTES9(A,B,Y,D,W2,NW2,INTERVAL) applies a 9 point Newton-Cotes rule to 
%   approximate the integral of F(X) from A to B.
%   The argument Y contains the 9 function values to be used in the
%   approximation and D is a constant to be used in the error estimate.
%   w2: is the rule's weights, nw2: is the nullrules weights
%   interval: logical signalling that too small intervals has been found.
%
%   See also COTES5
%   Modified by
%   T O Espelid, 10/07/2006

  h=(b-a)/2;x=linspace(a,b,9);Q=(h*w2)*y';
  Qabs=(h*abs(w2))*abs(y');
% define the local noise level
  noise=50*Qabs*eps;
%  Compute the error estimates
  e=(h*nw2)*y';e2=e.^2;E2= (e2(1:2:7)+e2(2:2:8));E=sqrt(E2);
  Emin=min(E(2:4));
  if Emin==0,
     rmax=2;
  else
     r=E(1:3)./E(2:4);
     if sum(isinf(r))>0
        rmax=2;
     else
        rmax=max(r);
     end
  end
  if rmax > 1
     err=D*max(E(1:4));
  elseif 0.25 < rmax
         err=D*rmax*E(1);
      else
         err=D*4*rmax*rmax*E(1);
  end
  
% If the highest degree null rules are on the local noise level, then put the error to zero.
  if (E(1)<noise) && (E(2) < noise), err=0; end;
%check if interval has become too small: unable to distinguish between the
%endpoints and the two points close to the endpoints.
  if ((x(2) <= a) || (b<=x(8)))
       interval =1;
%In order to avoid handling this interval again: put the error to zero.
       err=0;
   end

function [err,Q,Qabs]=cotes5(a,b,y,C,w1,nw1)
%
%
%   [err,Q,Qabs] = COTES5(A,B,Y,C,W1,NW1) applies a 5 point Newton-Cotes rule to 
%   approximate the integral of F(X) from A to B. w1 is the rule's weights
%   and nw1 the nullrules' weights. C is a constant. 
%   The argument Y contains the 5 function values to be used in the
%   approximation.
%
%   See also COTES9
%   Modified by
%   T O Espelid, 10/07/2006

  h=(b-a)/2;Q=(h*w1)*y';Qabs=(h*w1)*abs(y');
% define the local noise level 
  noise=50*Qabs*eps;E=abs((h*nw1)*y');
  Emin=min(E(2:4));
  if Emin==0,
     rmax=2;
  else
     r=E(1:3)./E(2:4);
     if sum(isinf(r))>0
        rmax=2;
     else
        rmax=max(r);
     end
  end
  if rmax > 1
     err=C*max(E(1:4));
  elseif 0.5 < rmax
         err=C*rmax*E(2);
      else
         err=8*C*(rmax^4)*E(2);
  end
% If the highest degree null rules are on the local noise level, then put the error to zero.
  if (E(1)<noise) && (E(2) < noise), err=0; end;

