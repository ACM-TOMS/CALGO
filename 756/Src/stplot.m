function [H,RE,IM] = stplot(w,beta,z,c,re,im,options)
%STPLOT Image of cartesian grid under Schwarz-Christoffel strip map.
%       STPLOT(W,BETA,Z,C) will adaptively plot the images under the
%       Schwarz-Christoffel exterior map of ten evenly spaced horizontal
%       and vertical lines in the upper half-plane. The abscissae of the
%       vertical lines will bracket the finite extremes of real(Z).  The
%       arguments are as in STPARAM.
%
%       STPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced
%       vertical and N evenly spaced horizontal lines.  Horizontal lines
%       are spaced to bracket real(Z); vertical lines are evenly spaced
%       between 0 and 1.
%
%       STPLOT(W,BETA,Z,C,RE,IM) will plot images of vertical lines
%       whose real parts are given in RE and horizontal lines whose
%       imaginary parts are given in IM.  Either argument may be empty.
%
%       STPLOT(W,BETA,Z,C,RE,IM,OPTIONS) allows customization of
%       HPPLOT's behavior.  See SCPLOTOPT.
%
%       H = STPLOT(W,BETA,Z,C,...) returns a vector of handles to all
%       the curves drawn in the interior of the polygon.  [H,RE,IM] =
%       STPLOT(W,BETA,Z,C,...) also returns the abscissae and ordinates
%       of the lines comprising the grid.
%	
%	See also SCPLOTOPT, STPARAM, STMAP, STDISP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

turn_off_hold = ~ishold;
N = length(w);
n = N-2;
w = w(:);
beta = beta(:);
z = z(:);
% Renumber vertices so that the ends of the strip map to w([1,k])
wend = [find(isinf(z)&(z<0)),find(isinf(z)&(z>0))];
renum = [wend(1):N,1:wend(1)-1];
w = w(renum);
beta = beta(renum);
z = z(renum);
k = find(renum==wend(2));
% nb = Number of prevertices on bottom edge of strip
nb = k-2;

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

minre = min(real(z(~isinf(z))));
maxre = max(real(z(~isinf(z))));
if (length(re)==1) & (re == round(re))
  if re < 1
    re = [];
  elseif re < 2
    re = mean([minre,maxre]);
  else
    m = re;
    re = linspace(minre,maxre,m);
    dre = diff(re(1:2));
    re = linspace(minre-dre,maxre+dre,m);
  end
end
if (length(im)==1) & (im == round(im))
  if im < 1
    im = [];
  else
    m = im;
    im = linspace(0,1,m+2);
    im([1,m+2]) = [];
  end
end

[nqpts,maxturn,maxlen,maxrefn] = scplotopt(options);

fig = gcf;
figure(fig);
plotpoly(w,beta);
drawnow
hold on

reflen = maxlen*max(abs(diff([w(~isinf(w));w(2)])));
qdat = scqdata(beta([2:k-1,k+1:N]),4);

for j = 1:length(re)
  zp = re(j) + i*linspace(0,1,15).';
  wp = stmap(zp,w,beta,z,c,qdat);
  bad = find(toobig(wp,maxturn,reflen,axis));
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    newz = [(zp(bad-1)+2*zp(bad))/3;(zp(bad+1)+2*zp(bad))/3];
    neww = stmap(newz,w,beta,z,c,qdat);
    [tmp,in] = sort(imag([zp;newz]));
    zp = [zp;newz];  wp = [wp;neww];
    zp = zp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig(wp,maxturn,reflen,axis));
  end
  linh(j) = plot(clipdata(wp,axis), 'g-','erasemode','none');
  drawnow
  set(linh(j),'erasemode','normal');
  Z(1:length(zp),j) = zp;
  W(1:length(wp),j) = wp;
end

x1 = min(-5,minre);
x2 = max(5,maxre);
axlim = axis;
for j = 1:length(im)
  zp = linspace(x1,x2,15).' + i*im(j);
  wp = stmap(zp,w,beta,z,c,qdat);
  bad = find(toobig([w(1);wp;w(k)],maxturn,reflen,axis))-1;
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    special = zeros(2,1);
    if isinf(w(1))
      if real(wp(1))>axlim(1) & real(wp(1))<axlim(2)&...
	    imag(wp(1))>axlim(3) & imag(wp(1))<axlim(4)
	special(1) = 1;
      end
    elseif any(bad==1)
      special(1) = 1;
    end
    if isinf(w(k))
      if real(wp(lenwp))>axlim(1) & real(wp(lenwp))<axlim(2)&...
	    imag(wp(lenwp))>axlim(3) & imag(wp(lenwp))<axlim(4)
	special(2) = 1;
      end
    elseif any(bad==lenwp)
      special(2) = 1;
    end
    bad(bad==1 | bad==lenwp) = [];
    newz = [(zp(bad-1)+2*zp(bad))/3;(zp(bad+1)+2*zp(bad))/3];
    zends = zp([1,lenwp]);
    newz = [newz;i*imag(zends(special))+5*real(zends(special))];
    neww = stmap(newz,w,beta,z,c,qdat);
    [tmp,in] = sort(real([zp;newz]));
    zp = [zp;newz];  wp = [wp;neww];
    zp = zp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig([w(1);wp;w(k)],maxturn,reflen,axis))-1;
  end
  linh(j+length(re)) = plot(clipdata([w(1);wp;w(k)],axis), ...
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

