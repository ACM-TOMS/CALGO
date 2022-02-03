function I = stquad(z1,z2,sing1,z,beta,qdat)
%STQUAD (not intended for calling directly by the user)
%	Numerical quadrature for the strip map.

%	z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%	integer indices which label the singularities in z1.  So if sing1(5)
%	= 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%	vector of *finite* singularities; beta is the vector of associated
%	turning angles.  z(1) must be the leftmost prevertex on the bottom
%	edge of the strip.  If nb=sum(~imag(z)), then z(1:nb) are on the
%	bottom edge and z(nb+1:n) are on the top (going right to left).
%	Note that length(beta) = length(z)+2, because the angles at the ends
%	of the strip are significant.  Hence beta(1) is the turn at the left
%	end and beta(nb+2) is at the right end.  qdat is quadrature data
%	from SCQDATA.
%
%	Make sure z and beta are column vectors.
%	
%	STQUAD integrates from a possible singularity at the left end to a
%	regular point at the right.  If both endpoints are singularities,
%	you must break the integral into two pieces and make two calls.
%	
%	The integral is subdivided, if necessary, so that no
%	singularity lies closer to the left endpoint than 1/2 the
%	length of the integration (sub)interval.
%
%	Written by Toby Driscoll.  Last updated 5/23/95.

n = length(z);
nb = sum(~imag(z));
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
  dist = min(1,2*min(abs(z([1:sng-1,sng+1:n])-za))/abs(zb-za));
  zr = za + dist*(zb-za);
  ind = rem(sng+n,n+1)+1;
  % Adjust Gauss-Jacobi nodes and weights to interval.
  nd = ((zr-za)*qdat(:,ind) + zr + za)/2; % G-J nodes
  wt = ((zr-za)/2) * qdat(:,ind+n+1); 	% G-J weights
  if any(~diff([za;nd;zr])) %| any(any(~terms)) 
    % Endpoints are practically coincident.
    I(k) = 0;
  else
    % Use Gauss-Jacobi on first subinterval, if necessary.
    if sng > 0
      wt = wt*(abs(zr-za)/2)^beta(sng+1+(sng>nb));
    end
    I(k) = stderiv(nd.',[-Inf;z(1:nb);Inf;z(nb+1:n)],beta,sng)*wt;
    while (dist < 1) & ~isnan(I(k))
      % Do regular Gaussian quad on other subintervals.
      zl = zr;
      dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
      zr = zl + dist*(zb-zl);
      nd = ((zr-zl)*qdat(:,n+1) + zr + zl)/2;
      wt = ((zr-zl)/2) * qdat(:,2*n+2);
      I(k) = I(k) + stderiv(nd.',[-Inf;z(1:nb);Inf;z(nb+1:n)],beta)*wt;
    end
  end
end

