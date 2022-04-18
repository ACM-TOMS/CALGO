function edgehandles = plotpoly(w,beta)
%PLOTPOLY Plot a (generalized) polygon.
%      PLOTPOLY(W,BETA) plots the polygon whose vertices are in vector W
%      and whose turning angles are in BETA.  Vertices at infinity are
%      permitted, but there must be at least two consecutive finite
%      vertices somewhere in W.
%
%	H = PLOTPOLY(W,BETA) returns a vector of handles to the polygon
%	sides.
%	
%	See also DRAWPOLY, MODPOLY.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

if nargin < 3
  thick = 3;
end
turn_off_hold = ~ishold;
%set(gcf,'defaultlinelinewidth',2.5*get(gcf,'defaultlinelinewidth'));
lw = 2*get(gcf,'defaultlinelinewidth');
n = length(w);
%w = w+i*1e-12;
wf = w(~isinf(w));
autoscale = strcmp(get(gca,'xlimmode'),'auto') & ...
    strcmp(get(gca,'ylimmode'),'auto');
autoscale = autoscale | turn_off_hold;
if autoscale
  lim = [min(real(wf)),max(real(wf)),min(imag(wf)),max(imag(wf))];
  maxdiff = max(diff(lim(1:2)),diff(lim(3:4)));
else
  lim = axis;
end
first = 1;
if ~any(isinf(w))
  for j = 1:n-1
    edgeh(j) = plot(real(w(j:j+1)),imag(w(j:j+1)),'-','linewid',lw);
    if j==1, hold on; end;
  end
  edgeh(n) = plot(real(w([n,1])), imag(w([n,1])),'-','linewid',lw);
  if autoscale
    lim(1:2) = mean(lim(1:2)) + 0.55*maxdiff*[-1,1];
    lim(3:4) = mean(lim(3:4)) + 0.55*maxdiff*[-1,1];
  end
  axis(lim)
else
  if any(isinf(w(1:2)))
    first = min(find(~isinf(w) & ~isinf(w([2:n,1]))));
    if isempty(first), 
      error('There must be two consecutive finite vertices.')
    end
    w = w([first:n,1:first-1]);
    beta = beta([first:n,1:first-1]);
  end
  edgeh(1) = plot(real(w(1:2)),imag(w(1:2)),'-','linewid',lw);
  ang = angle(w(2)-w(1));
  if autoscale
    lim(1:2) = mean(lim(1:2)) + 0.65*maxdiff*[-1,1];
    lim(3:4) = mean(lim(3:4)) + 0.65*maxdiff*[-1,1];
  end
  R = max(lim(2)-lim(1),lim(4)-lim(3));
  axis(lim)
  hold on
  j = 2;
  while j < n
    if ~isinf(w(j+1))
      edgeh(j) = plot(real(w(j:j+1)),imag(w(j:j+1)),'-','linewid',lw);
      ang = ang - pi*beta(j);
      j = j+1;
    else
      ang = ang-pi*beta(j);
      z = [w(j);w(j)+R*exp(i*ang)];
      edgeh(j) = plot(real(z),imag(z),'-','linewid',lw);
      ang = ang-pi*beta(j+1);
      z = [w(rem(j+1,n)+1)-R*exp(i*ang);w(rem(j+1,n)+1)];
      edgeh(j+1) = plot(real(z),imag(z),'-','linewid',lw);
      j = j+2;
    end
  end 
  if j==n
    edgeh(n) = plot(real(w([n,1])),imag(w([n,1])),'-','linewid',lw);
  end
end

axis square
axis equal
if nargout 
  edgehandles([first:n,1:first-1]) = edgeh;
end
%set(gcf,'defaultlinelinewidth',get(gcf,'defaultlinelinewidth')/1.6);
if turn_off_hold
  hold off
end 

