function [H,R2,THETA] = deplot(w,beta,z,c,R,theta,options)
%DEPLOT Image of polar grid under Schwarz-Christoffel exterior map.
%       DEPLOT(W,BETA,Z,C) will adaptively plot the images under the
%       Schwarz-Christoffel exterior map of ten evenly spaced circles
%       and rays in the unit disk.  The arguments are as in DEPARAM.
%
%       DEPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced
%       circles and N evenly spaced rays.
%
%       DEPLOT(W,BETA,Z,C,R,THETA) will plot images of circles whose
%       radii are given in R and rays whose arguments are given in
%       THETA.  Either argument may be empty.
%
%       DEPLOT(W,BETA,Z,C,R,THETA,OPTIONS) allows customization of
%       DEPLOT's behavior.  See SCPLOTOPT.
%
%       H = DEPLOT(W,BETA,Z,C,...) returns a vector of handles to all
%       the curves drawn in the interior of the polygon.  [H,R,THETA] =
%       DEPLOT(W,BETA,Z,C,...) also returns the moduli and arguments of
%       the curves comprising the grid.
%	
%	See also SCPLOTOPT, DEPARAM, DEMAP, DEDISP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
beta = beta(:);
z = z(:);

turn_off_hold = ~ishold;
if nargin < 7
  options = [];
  if nargin < 6
    theta = [];
    if nargin < 5 
      R = [];
    end
  end
end

if isempty([R(:);theta(:)])
  R = 10;
  theta = 10;
end

if (length(R)==1) & (R == round(R))
  m = R+2;
  R = fliplr(linspace(.25,1,m));
  R([1,m]) = [];
end
if (length(theta)==1) & (theta == round(theta))
  m = theta+1;
  theta = linspace(0,2*pi,m);
  theta(m) = [];
end

[nqpts,maxturn,maxlen,maxrefn] = scplotopt(options);
autoscale = strcmp(get(gca,'xlimmode'),'auto') & ...
    strcmp(get(gca,'ylimmode'),'auto');
autoscale = autoscale | ~ishold;

fig = gcf;
figure(fig);
plotpoly(w,beta);
hold on

axlim = axis;
if autoscale
  axlim(1:2) = axlim(1:2) + 0.25*diff(axlim(1:2))*[-1,1];
  axlim(3:4) = axlim(3:4) + 0.25*diff(axlim(3:4))*[-1,1];
  axis(axlim);
end
drawnow

n = length(w);
wf = w(~isinf(w));
reflen = maxlen*max(abs(diff([wf;wf(1)])));

qdat = scqdata(beta,nqpts);
Rp0 = linspace(.1,1,15)';
  
for j = 1:length(R)
  tp = linspace(0,2*pi,16)';
  tp = [tp(length(tp)-1)-2*pi;tp];
  zp = R(j)*exp(i*tp);
  wp = demap(zp,w,beta,z,c,qdat);
  bad = find(toobig(wp,maxturn,reflen,axis));
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    newt = [(tp(bad-1)+2*tp(bad))/3;(tp(bad+1)+2*tp(bad))/3];
    newz = R(j)*exp(i*newt);
    neww = demap(newz,w,beta,z,c,qdat);
    [k,in] = sort([tp;newt]);
    tp = [tp;newt];  wp = [wp;neww];
    tp = tp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig(wp,maxturn,reflen,axis));
  end
  tp(tp<0) = tp(tp<0) + 2*pi;
  [k,in] = sort(tp);
  linh(j) = plot(clipdata(wp(in),axis), 'g-','erasemode','none');
  set(linh(j),'erasemode','normal');
  drawnow
  Z(1:length(zp),j) = zp;
  W(1:length(wp),j) = wp;
end

for j = 1:length(theta)
  Rp = Rp0;
  zp = Rp*exp(i*theta(j));
  wp = demap(zp,w,beta,z,c,qdat);
  bad = find(toobig(wp,maxturn,reflen,axis));
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    newR = [(Rp(bad-1)+2*Rp(bad))/3;(Rp(bad+1)+2*Rp(bad))/3];
    newz = newR*exp(i*theta(j));
    neww = demap(newz,w,beta,z,c,qdat);
    [k,in] = sort([Rp;newR]);
    Rp = [Rp;newR];  wp = [wp;neww];
    Rp = Rp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig(wp,maxturn,reflen,axis));
  end
  linh(j+length(R)) = plot(clipdata(wp,axis), 'g-','erasemode','none');
  drawnow
  set(linh(j+length(R)),'erasemode','normal');
  Z(1:length(zp),j+length(R)) = zp;
  W(1:length(wp),j+length(R)) = wp;
end

% Force redraw to get clipping enforced.
set(fig,'color',get(fig,'color'))
if turn_off_hold, hold off, end;
if nargout > 0
  H = linh;
  if nargout > 1
    R2 = R;
    if nargout > 2
      THETA = theta;
    end
  end
end 

