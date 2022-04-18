function bad = toobig(wp,maxturn,reflen,clip)
%TOOBIG (not intended for calling directly by the user)
%	Used by plotting functions to spot points in a polyline at which
%	turning angles or line segments are too large.
%
%       See also HPPLOT, DPLOT, DEPLOT, STPLOT, RPLOT.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

m = length(wp);
dwp = abs(diff(wp(:)));
dwp = max([dwp(1:m-2).';dwp(2:m-1).']).';
angl = abs(scangle(wp(:)));
bad = [0;(((angl(2:m-1) > maxturn/180)&(dwp > reflen/8)) | (dwp > reflen));0];
if nargin >= 4				% ignore clipped points
  xp = real(wp);
  yp = imag(wp);
  inside = (xp>clip(1)) & (xp<clip(2)) & (yp>clip(3)) & (yp<clip(4));
  bad = bad & (inside | [1;inside(1:m-1)] | [inside(2:m);1]);
end
