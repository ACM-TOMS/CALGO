% STIELTJES Discretized Stieltjes procedure.
%
%    Given the discrete inner product whose nodes are contained 
%    in the first column, and whose weights are contained in the
%    second column, of the nx2 array xw, the call ab=STIELTJES(n,xw) 
%    generates the first n recurrence coefficients ab of the 
%    corresponding discrete orthogonal polynomials. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab.
%
function ab=stieltjes(N,xw)
tiny=10*realmin;
huge=.1*realmax;
%
% Remove data with zero weights
%
xw=sortrows(xw,2);
index=min(find(xw(:,2)~=0));
xw=xw(index:size(xw(:,2)),:);
xw=sortrows(xw,1);
Ncap=size(xw,1);
%
if(N<=0|N>Ncap), error('N in sti out of range'), end
s0=ones(1,Ncap)*xw(:,2);
ab(1,1)=xw(:,1)'*xw(:,2)/s0; ab(1,2)=s0;
if N==1, return, end
p1=zeros(Ncap,1); p2=ones(Ncap,1);
%
% The following scaling has to be adopted in case the
% orthogonal polynomials underflow or overflow. In the
% former case, c has to be chosen sufficiently large 
% positive, in the latter case sufficiently small
% positive. (The example below relates to Table2_9.m.)
%
%if N==320
%  c=1e50;
%  p2=c*p2;
%  s0=c^2*s0;
%end
%
for k=1:N-1
  p0=p1; p1=p2;
  p2=(xw(:,1)-ab(k,1)).*p1-ab(k,2)*p0;
  s1=xw(:,2)'*(p2.^2);
  s2=xw(:,1)'*(xw(:,2).*(p2.^2));
  if(max(abs(p2))>huge)|(abs(s2)>huge)
    error(sprintf('impending overflow in stieltjes for k=%3.0f',k))
  end
  if abs(s1)<tiny
    error(sprintf('impending underflow in stieltjes for k=%3.0f',k))
  end
  ab(k+1,1)=s2/s1; ab(k+1,2)=s1/s0;
  s0=s1;
end
