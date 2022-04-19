function I = hpquad(z1,z2,sing1,x,beta,qdat)
%HPQUAD (not intended for calling directly by the user)
%	Numerical quadrature for the half-plane map.

%	z1,z2 are vectors of left and right endpoints.  sing1 is a vector
%	of integer indices which label the singularities in z1.  So if
%	sing1(5) = 3, then z1(5) = x(3).  A zero means no singularity. 
%	x is the vector of finite singularities;  beta is the vector of
%	associated turning angles.  qdat is quadrature data from SCQDATA.
%
%	Make sure x and beta are column vectors.
%	
%	HPQUAD integrates from a possible singularity at the left end to a
%	regular point at the right.  If both endpoints are singularities,
%	you must break the integral into two pieces and make two calls.
%	
%	The integral is subdivided, if necessary, so that no
%	singularity lies closer to the left endpoint than 1/2 the
%	length of the integration (sub)interval.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

nqpts = size(qdat,1);
% Note: Here n is the total number of *finite* singularities; i.e., the
% number of terms in the product appearing in the integrand.
n = length(x);
bigx = x(:,ones(1,nqpts));
bigbeta = beta(:,ones(1,nqpts));
if isempty(sing1)
  sing1 = zeros(length(z1),1);
end

I = zeros(size(z1));
nontriv = find(z1(:)~=z2(:))';

for k = nontriv
  za = z1(k);
  zb = z2(k);
  sng = sing1(k);

  % Allowable integration step, based on nearest singularity.
  dist = min(1,2*min(abs(x([1:sng-1,sng+1:n])-za))/abs(zb-za));
  zr = za + dist*(zb-za);
  ind = rem(sng+n,n+1)+1;
  % Adjust Gauss-Jacobi nodes and weights to interval.
  nd = ((zr-za)*qdat(:,ind) + zr + za)/2; % G-J nodes
  wt = ((zr-za)/2) * qdat(:,ind+n+1); 	% G-J weights
  terms = nd(:,ones(n,1)).' - bigx;
  if any(~diff(nd)) | any(any(~terms)) 
    % Endpoints are practically coincident.
    I(k) = 0;
  else
    % Use Gauss-Jacobi on first subinterval, if necessary.
    if sng > 0
      terms(sng,:) = terms(sng,:)./abs(terms(sng,:));
      wt = wt*(abs(zr-za)/2)^beta(sng);
    end
    I(k) = exp(sum(log(terms).*bigbeta))*wt;
    while dist < 1              
      % Do regular Gaussian quad on other subintervals.
      zl = zr;
      dist = min(1,2*min(abs(x-zl))/abs(zl-zb));
      zr = zl + dist*(zb-zl);
      nd = ((zr-zl)*qdat(:,n+1) + zr + zl)/2;
      wt = ((zr-zl)/2) * qdat(:,2*n+2);
      terms = nd(:,ones(n,1)).' - bigx;
      I(k) = I(k) + exp(sum(log(terms).*bigbeta)) * wt;
    end
  end
end

