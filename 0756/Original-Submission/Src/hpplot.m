function [H,RE,IM] = hpplot(w,beta,x,c,re,im,options)
%HPPLOT Image of cartesian grid under Schwarz-Christoffel half-plane map.
%       HPPLOT(W,BETA,X,C) will adaptively plot the images under the
%       Schwarz-Christoffel exterior map of ten evenly spaced horizontal
%       and vertical lines in the upper half-plane. The abscissae of the
%       vertical lines will bracket the finite extremes of X.  The
%       arguments are as in HPPARAM.
%
%       HPPLOT(W,BETA,X,C,M,N) will plot images of M evenly spaced
%       vertical and N evenly spaced horizontal lines.  The spacing will
%       be the same in both directions.
%
%       HPPLOT(W,BETA,X,C,RE,IM) will plot images of vertical lines
%       whose real parts are given in RE and horizontal lines whose
%       imaginary parts are given in IM.  Either argument may be empty.
%
%       HPPLOT(W,BETA,X,C,RE,IM,OPTIONS) allows customization of
%       HPPLOT's behavior.  See SCPLOTOPT.
%
%       H = HPPLOT(W,BETA,X,C,...) returns a vector of handles to all
%       the curves drawn in the interior of the polygon.  [H,RE,IM] =
%       HPPLOT(W,BETA,X,C,...) also returns the abscissae and ordinates
%       of the lines comprising the grid.
%	
%	See also SCPLOTOPT, HPPARAM, HPMAP, HPDISP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

turn_off_hold = ~ishold;
n = length(w);
w = w(:);
beta = beta(:);
x = x(:);
if nargin < 7
  options = [];
  if nargin < 6
    im = [];
    if nargin < 5
      re = [];
    end
  end
end

if isempty([re(:);im(:)])
  re = 10;
  im = 10;
end

if (length(re)==1) & (re == round(re))
  if re < 1
    re = [];
  elseif re < 2
    re = mean(x([1,n-1]));
  else
    m = re;
    re = linspace(x(1),x(n-1),m);
    dre = diff(re(1:2));
    re = linspace(x(1)-dre,x(n-1)+dre,m);
  end
end
if (length(im)==1) & (im == round(im))
  if length(re) < 2
    im = linspace(0,4,im+1);
    im(1) = [];
  else
    im = mean(diff(re))*(1:im);
  end
end

[nqpts,maxturn,maxlen,maxrefn] = scplotopt(options);

fig = gcf;
figure(fig);
plotpoly(w,beta);
drawnow
hold on

n = length(w);
reflen = maxlen*max(abs(diff([w(~isinf(w));w(1)])));
if any(isinf(x))
  qdat = scqdata(beta(1:n-1),nqpts);
else
  qdat = scqdata(beta,nqpts);
end

y2 = max(x(n-1),10);
for j = 1:length(re)
  zp = re(j) + i*linspace(0,y2,15).';
  wp = hpmap(zp,w,beta,x,c,qdat);
  bad = find(toobig([wp;w(n)],maxturn,reflen,axis));
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    newz = [];
    special = find(bad==lenwp);
    newz = re(j) + i*5*imag(zp(bad(special)));
    bad(special) = [];
    newz = [newz;(zp(bad-1)+2*zp(bad))/3;(zp(bad+1)+2*zp(bad))/3];
    neww = hpmap(newz,w,beta,x,c,qdat);
    [k,in] = sort(imag([zp;newz]));
    zp = [zp;newz];  wp = [wp;neww];
    zp = zp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig([wp;w(n)],maxturn,reflen,axis));
  end
  linh(j) = plot(clipdata([wp;w(n)],axis), 'g-','erasemode','none');
  drawnow
  set(linh(j),'erasemode','normal');
  Z(1:length(zp),j) = zp;
  W(1:length(wp),j) = wp;
end

x1 = min(-10,x(n-1));
x2 = max(40,x(n-1));
axlim = axis;
for j = 1:length(im)
  zp = linspace(x1,x2,15).' + i*im(j);
  wp = hpmap(zp,w,beta,x,c,qdat);
  bad = find(toobig([w(n);wp;w(n)],maxturn,reflen,axis)) - 1;
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    special = zeros(2,1);
    if isinf(w(n))
      ends = wp([1,lenwp]);
      special = real(ends)>axlim(1) & real(ends)<axlim(2) & ...
	  imag(ends)>axlim(3) & imag(ends)<axlim(4);
    else 
      special(1) = any(bad==1);
      special(2) = any(bad==lenwp);
    end
    bad(bad==1 | bad==lenwp) = [];
    newz = [(zp(bad-1)+2*zp(bad))/3;(zp(bad+1)+2*zp(bad))/3];
    zends = zp([1,lenwp]);
    newz = [newz;i*imag(zends(special))+5*real(zends(special))];
    neww = hpmap(newz,w,beta,x,c,qdat);
    [k,in] = sort(real([zp;newz]));
    zp = [zp;newz];  wp = [wp;neww];
    zp = zp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig([w(n);wp;w(n)],maxturn,reflen,axis)) - 1;
  end
  linh(j+length(re)) = plot(clipdata([w(n);wp;w(n)],axis),...
      'g-','erasemode','none');
  drawnow
  set(linh(j+length(re)),'erasemode','normal');
  Z(1:length(zp),j+length(re)) = zp;
  W(1:length(wp),j+length(re)) = wp;
end
  
% Force redraw to get clipping enforced.
set(fig,'color',get(fig,'color'))
if turn_off_hold, hold off, end;
if nargout > 0
  H = linh;
  if nargout > 1
    RE = re;
    if nargout > 2
      IM = im;
    end
  end
end 

