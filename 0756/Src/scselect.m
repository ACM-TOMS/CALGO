function K = scselect(w,beta,m)
%SCSELECT Select one or more vertices in a polygon.
%	K = SCSELECT(W,BETA,M) draws the polygon given by W and BETA
%	into the current figure window and then allows the user to
%	select M vertices using the mouse.  If M is not given, it
%	defaults to 1.  On exit K is a vector of indices into W.
%
%       See also DRAWPOLY, PLOTPOLY, MODPOLY.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
if any(isinf(w) & isinf(w([2:n,1])))
  error('Infinite vertices must not be adjacent')
end

plotpoly(w,beta)
turn_off_hold = ~ishold;
hold on

first = min(find(~isinf(w) & ~isinf(w([2:n,1]))));
renum = [first:n,1:first-1];
w = w(renum);
beta = beta(renum);

axlim = axis;
maxdiff = max(diff(axlim(1:2)),diff(axlim(3:4)));
axlim(1:2) = mean(axlim(1:2)) + 0.57*maxdiff*[-1.05,1];
axlim(3:4) = mean(axlim(3:4)) + 0.57*maxdiff*[-1,1];
    
h = zeros(n,2);
h(1,1) = plot(real(w(1)),imag(w(1)),'.','mark',22);
ang = angle(w(2)-w(1));
colrs = get(gca,'colororder');
colr = colrs(1,:);
for j = 2:n
  if ~isinf(w(j))
    if ~imag(w(j))
      w(j) = w(j) + eps*i;
    end 
    h(j,1) = plot(w(j),'.','mark',22);
    ang = ang - pi*beta(j);
  else
    for p = 1:2
      theta = ang + pi*(p==2);
      base = w(rem(j-2+2*(p==2),n)+1);
      Rx = (axlim(1:2) - real(base)) / (cos(theta)+eps*(cos(theta)==0));
      Ry = (axlim(3:4) - imag(base)) / (sin(theta)+eps*(sin(theta)==0));
      R = [Rx,Ry];
      wj = base + min(R(R>0))*exp(i*theta);
      str = sprintf('inf (%i)',renum(j));
      h(j,p) = text(real(wj),imag(wj),str,'horiz','center',...
	  'fontsize',14,'fontweight','bold','color',colr);
      ang = ang - pi*beta(j)*(p==1);
    end 
  end 
end 

colr = colrs(min(2,size(colrs,1)),:);
oldptr = get(gcf,'pointer');
set(gcf,'pointer','circle');
if nargin < 3
  m = 1;
end
for j = 1:m
  k = [];
  while isempty(k)
    waitforbuttonpress;
    obj = get(gcf,'currentobj');
    [k,tmp] = find(obj==h);
    if isempty(k)
      disp('Selected object not a vertex.  Try again.')
    end
  end
  set(h(k,(h(k,:)>0)),'color',colr)
  drawnow
  K(j) = k;
end
set(gcf,'pointer',oldptr)

delete(h(h>0))
drawnow
if turn_off_hold
  hold off
end 

K = renum(K);
