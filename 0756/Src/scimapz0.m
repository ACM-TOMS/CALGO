function [z0,w0] = scimapz0(prefix,wp,w,beta,z,c,qdat,aux)
%SCIMAPZ0 (not intended for calling directly by the user)
%       SCIMAPZ0 returns starting points for computing inverses of
%       Schwarz-Christoffel maps.
%
%	Each wp(j) (in the polygon plane) requires z0(j) (in the
%	fundamental domain) whose image w0(j) is such that the line
%	segment from w0(j) to wp(j) lies in the target (interior or
%	exterior) region.  The algorithm here is to choose z0(j) as a
%	(weighted) average of successive pairs of adjacent prevertices.
%	The resulting w0(j) is on a polygon side.  Each choice is tested
%	by looking for intersections of the segment with (other) sides
%	of the polygon.
%
%	After randomly trying 10 weights with such prevertex pairs, the
%	routine gives up.  Failures are pretty rare.  Slits are the most
%	likely cause of trouble, since the intersection method doesn't
%	know "which side" of the slit it's on.  In such a case you will
%	have to supply starting points manually, perhaps by a
%	continuation method.
%
%	See also HPINVMAP, DINVMAP, DEINVMAP, RINVMAP, STINVMAP.
%
%	Written by Toby Driscoll.  Last updated 7/7/95.

%	P.S. This file illustrates why the different domains in the SC
%	Toolbox have mostly independent M-files.  The contingencies for
%	the various geometries become rather cumbersome.

n = length(w);
shape = wp;
wp = wp(:);
z0 = wp;
w0 = wp;
from_disk = strcmp(prefix(1),'d');
from_hp = strcmp(prefix,'hp');
from_strip = strcmp(prefix,'st');
from_rect = strcmp(prefix,'r');
if from_strip
  kinf = max(find(isinf(z)));
  argw = cumsum([angle(w(3)-w(2));-pi*beta([3:n,1])]);
  argw = argw([n,1:n-1]);
else
  argw = cumsum([angle(w(2)-w(1));-pi*beta(2:n)]);
end
if from_disk
  argz = angle(z);
  argz(argz<=0) = argz(argz<=0) + 2*pi;
end

factor = 0.5;				% images of midpoints of preverts
done = zeros(1,length(wp));
m = length(wp);
iter = 0;

while m > 0				% while some not done
  % Choose a point on each side of the polygon.
  for j = 1:n 	
    if from_disk
      if j<n
	zbase(j) = exp(i*(factor*argz(j) + (1-factor)*argz(j+1)));
      else 
	zbase(j) = exp(i*(factor*argz(n) + (1-factor)*(2*pi+argz(1))));
      end
    elseif from_hp
      if j < n-1			% between two finite points
	zbase(j) = z(j) + factor*(z(j+1)-z(j));
      elseif j==n-1			% between x(n-1) & Inf
	zbase(j) = max(10,z(n-1))/factor;
      else				% between -Inf and x(1)
	zbase(j) = min(-10,z(1))/factor;
      end
    elseif from_strip
      if j==1
	zbase(j) = min(-1,real(z(2)))/factor/4;
      elseif j==kinf-1
	zbase(j) = max(1,real(z(kinf-1)))/factor/4;
      elseif j==kinf
	zbase(j) = i+max(1,real(z(kinf+1)))/factor/4;
      elseif j==n
	zbase(j) = i+min(-1,real(z(n)))/factor/4;
      else 
	zbase(j) = z(j) + factor*(z(j+1)-z(j));
      end
    elseif from_rect
      zbase(j) = z(j) + factor*(z(rem(j,n)+1)-z(j));
    end
    if ~from_rect
      wbase(j) = feval([prefix,'map'],zbase(j),w,beta,z,c,qdat);
    else 
      wbase(j) = feval([prefix,'map'],zbase(j),w,beta,z,c,qdat,aux);
    end
    
  end

  % Now, cycle thru starting points
  for j = 1:n 				
    z0(~done) = ones(m,1)*zbase(j);
    w0(~done) = ones(m,1)*wbase(j);
    notdone = find(~done);
    done = ones(1,length(wp));
    % Test line segment for intersections with other sides.
    % We'll parameterize line segment and polygon side, compute parameters
    % at intersection, and check parameters at intersection.
    for k=[1:j-1,j+1:n]
      if isinf(w(k))
	A(:,1) = [cos(argw(k)+pi);sin(argw(k)+pi)];
	wk = w(rem(k,n)+1);
	s1max = Inf;
      else
	A(:,1) = [cos(argw(k));sin(argw(k))];
	wk = w(k);
	s1max = abs(w(rem(k,n)+1)-w(k));
      end
      for p = notdone
	A(:,2) = [real(w0(p)-wp(p));imag(w0(p)-wp(p))];
	% Get line segment and side parameters at intersection.
	s = A\[real(w0(p)-wk);imag(w0(p)-wk)];
	% Intersection occurs interior to side? and segment? 
	if s(1)>=0 & s(1)<=s1max 
	  if abs(s(2)-1) < 30*eps
	    % Special case: wp(p) is on polygon side k
	    z0(p) = zbase(k);
	    w0(p) = wbase(k);
	  elseif s(2) > -10*eps & s(2) < 1
	    % Intersection interior to segment: it's no good
	    done(p) = 0;
	  end
	end
      end
    end
    m = sum(~done);
    if ~m, break, end
  end
  if iter > 10
    error('Can''t seem to choose starting points.  Supply them yourself.')
  else
    iter = iter + 1;
  end
  factor = rand(1);			% abandon midpoints
end

shape(:) = z0;
z0 = shape;
shape(:) = w0;
w0 = shape;
