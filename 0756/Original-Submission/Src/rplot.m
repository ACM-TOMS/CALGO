function [H,RE,IM] = rplot(w,beta,z,c,L,re,im,options)
%RPLOT  Image of cartesian grid under Schwarz-Christoffel rectangle map.
%       RPLOT(W,BETA,Z,C,L) will adaptively plot the images under the
%       Schwarz-Christoffel rectangle map of ten evenly spaced
%       horizontal and vertical lines in the retangle RECT.  The
%       arguments are as in RPARAM.
%
%       RPLOT(W,BETA,Z,C,L,M,N) will plot images of M evenly spaced
%       vertical and N evenly spaced horizontal lines.
%
%       RPLOT(W,BETA,Z,C,L,RE,IM) will plot images of vertical lines
%       whose real parts are given in RE and horizontal lines whose
%       imaginary parts are given in IM.  Either argument may be empty.
%
%       RPLOT(W,BETA,Z,C,L,RE,IM,OPTIONS) allows customization of
%       RPLOT's behavior.  See SCPLOTOPT.
%
%       H = RPLOT(W,BETA,Z,C,L,...) returns a vector of handles to all
%       the curves drawn in the interior of the polygon.  [H,RE,IM] =
%       RPLOT(W,BETA,Z,C,L,...) also returns the abscissae and ordinates
%       of the lines comprising the grid.
%	
%	See also SCPLOTOPT, RPARAM, RMAP, RDISP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

turn_off_hold = ~ishold;
n = length(w);
w = w(:);
beta = beta(:);
z = z(:);
[w,beta,z,corners] = rcorners(w,beta,z);
rect = z(corners);

if nargin < 8
  options = [];
  if nargin < 7
    im = [];
    if nargin < 7
      re = [];
    end
  end
end

Kp = imag(rect(2));
K = rect(1);

if isempty([re(:);im(:)])
  re = 10;
  im = 10;
end

if (length(re)==1) & (re == round(re))
  if re < 1
    re = [];
  else
    m = re;
    re = linspace(-K,K,m+2);
    re([1,m+2]) = [];
  end
end
if (length(im)==1) & (im == round(im))
  if im < 1
    im = [];
  else
    m = im;
    im = linspace(0,Kp,m+2);
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
qdat = scqdata(beta,4);

for j = 1:length(re)
  zp = re(j) + i*linspace(0,Kp,15).';
  wp = rmap(zp,w,beta,z,c,L,qdat);
  bad = find(toobig(wp,maxturn,reflen,axis));
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    newz = [(zp(bad-1)+2*zp(bad))/3;(zp(bad+1)+2*zp(bad))/3];
    neww = rmap(newz,w,beta,z,c,L,qdat);
    [k,in] = sort(imag([zp;newz]));
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

for j = 1:length(im)
  zp = linspace(-K,K,15).' + i*im(j);
  wp = rmap(zp,w,beta,z,c,L,qdat);
  bad = find(toobig(wp,maxturn,reflen,axis));
  iter = 0;
  while (~isempty(bad)) & (iter < maxrefn)
    lenwp = length(wp);
    newz = [(zp(bad-1)+2*zp(bad))/3;(zp(bad+1)+2*zp(bad))/3];
    neww = rmap(newz,w,beta,z,c,L,qdat);
    [k,in] = sort(real([zp;newz]));
    zp = [zp;newz];  wp = [wp;neww];
    zp = zp(in);     wp = wp(in);
    iter = iter + 1;
    bad = find(toobig(wp,maxturn,reflen,axis));
  end
  linh(j+length(re)) = plot(clipdata(wp,axis), 'g-','erasemode','none');
  drawnow
  Z(1:length(zp),j+length(re)) = zp;
  W(1:length(wp),j+length(re)) = wp;
  set(linh(j+length(re)),'erasemode','normal');
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

