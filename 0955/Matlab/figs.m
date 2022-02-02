%
% first figure showing various regions in (x,lambda) space
%

function figs()

close all
clear all

figure(1)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[0.8 0.8]; set(gcf,'pos',pos);

x1 = 1;
x2 = 10;
x3 = 1000;
x4 = 2.5;
x4 = 4;

%
% bottom-up region
%

loglog([x2 x2],[x4 x3],'-.k',...
       [x2 x3],[x4 x4],'-.k',...
       [x1 x2],[x1 x2],'-.k'); hold on
axis square; axis([x1 x3 x1 x3]);
xlabel('x'); ylabel('$\lambda$','interpreter','latex'); 
text(1.5,50,'bottom-up')
text(1.5,1.25,'bottom-up')
text(2.0,1.75,'top-down')
text(12,16,'|w|<3')

%
% r=x/lambda lines for central Temme region
%

rhi = 4.0;
rhi = 3.25;
rlo = 0.4;
loglog([x4*rhi x3],[x4 x3/rhi],'--k', ...
       [x2 x3*rlo],[x2/rlo x3],'--k');

%
% curves corresponding to different values of u
%

x = exp(linspace(log(x1),log(x3),1000));

%for u = [3e-4 -3e-4 6e-8 -6e-8 1e-300 -eps(0.5)]
%  for u = [ 3e-4   -3e-4    1e-300 -eps(0.5)  1e-38  -2^(-24) ]

% for u = [ 3e-4 -3e-4 1e-38 -2^(-24) 1e-300 -eps(0.5) ]

for u = [ 3e-4 -3e-4 1e-38 -2^(-24) 1e-300 -1e-300 ]

  if (u>0)
    y = gammaincinv(u,x,'upper');
  else
%    y = gammaincinv(-u,x,'lower');

    for i = 1:length(x)
      xx = x(i);

      yy  = xx;
      fac = 2;
      for kount = 1:20
         while(gammainc(yy/fac,xx,'lower')>-u)
           yy = yy/fac;
         end
         fac = sqrt(fac);
      end

      y(i) = yy;
    end

  end
  loglog(x,y,'k');
end

print('-deps2','fig1.eps')

% return

%
% second figure illustrating CDF inversion
%

figure(2)
%pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.8]; set(gcf,'pos',pos);
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.2 0.8]; set(gcf,'pos',pos);

N  = 18;
lam = 10;

x = linspace(4,N);
plot(x,gammainc(lam,x,'upper'),'k'); hold on

n = 4:N;
n = 0:N;
p = poisscdf(n,lam);

ns = reshape([n; n],2*length(n),1);
ps = reshape([p; p],2*length(n),1);
plot(ns(2:end),ps(1:end-1),'k'); hold on
xlabel('x'); ylabel('u')

x  = 10.8;
u  = gammainc(lam,x,'upper');
x2 = poissinv(u,lam);

[ax, ay] = dsxy2figxy(gca, [n(1) x2],[u u]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [x2 x2],[u 0]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [n(1) x],[u u]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [x x],[u 0]);
annotation('arrow',ax,ay,'HeadWidth',6)
[ax, ay] = dsxy2figxy(gca, [x x2],[0 0]);
annotation('arrow',ax,ay,'HeadWidth',6)

print('-deps2','fig2.eps')


%
% 4 figures illustrating accuracy of Normal approximations
% -- last one not used in paper
%

figure(3)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.5]; set(gcf,'pos',pos);
figure(4)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.5]; set(gcf,'pos',pos);
figure(5)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.0]; set(gcf,'pos',pos);

xs  = [10 100 exp(linspace(log(10),log(1000),1000))];

for k = 1:length(xs)
  x = xs(k);
  w = linspace(-3,3,1000);

  lam = gammaincinv(normcdf(w),x,'upper');

  x1 = lam + sqrt(lam).*w;
  x2 = x1 + w.^2/6 + 1/3;
  x3 = x2 + (- w.^3/72 - w/36)./sqrt(lam);
  x4 = x3 + (w.^4/270 + (7*w.^2)/810 - 8/405)./lam;
  delta =   (w.^4/160 +    w.^2 /80  + 1/40 )./lam;

  if k<3
    figure(3)
    subplot(3,2,k); plot(lam,x2-x,'k')
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 1')
    title(['x = ' num2str(x)])
    subplot(3,2,k+2); plot(lam,x3-x,'k')
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 2')
    subplot(3,2,k+4); plot(lam,x4-x,'k')
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 3')

    figure(4)
    subplot(1,2,k);   plot(lam,(x3-x)./delta,'k')
    xlabel('$\lambda$','interpreter','latex'); ylabel('relative error 2')
    title(['x = ' num2str(x)])
  else
    e1(k-2) = max(abs(x1-x));
    e2(k-2) = max(abs(x2-x));
    e3(k-2) = max(abs(x3-x));
    e4(k-2) = max(abs(x4-x));
    re3(k-2) = max(abs(x3-x)./delta);
  end
end

figure(3)
print('-deps2','asymp1.eps');

figure(4)
print('-deps2','asymp3.eps');

figure(5)
loglog(xs(3:end),e2,'-.k',xs(3:end),e3,'--k',xs(3:end),e4,'--k')
axis([10 1000 1e-4 1]) 
xlabel('x')
ylabel('maximum error')
legend('error 1','error 2','error 3','Location','Southwest')
print('-deps2','asymp2.eps');

figure(6)
semilogx(xs(3:end),re3,'-k')
axis([10 1000 0 1]) 
xlabel('x'); ylabel('maximum relative error')


%
% final figure plotting accuracy results computed by poissinv_check
%

res = [ ...
         1    8.94e-08    5.55e-17  5.55e-17  2.02e-16 ; ...
         2    8.94e-08    1.94e-16  1.94e-16  3.06e-16 ; ...
         4    1.21e-07    1.11e-16  1.70e-16  3.13e-16 ; ...
         8    5.19e-07    3.86e-16  3.37e-16  4.56e-16 ; ...
        16    7.05e-07    7.41e-16  9.51e-16  7.45e-16; ...
        32    7.61e-07    1.23e-15  1.82e-15  1.87e-15 ; ...
        64    6.58e-07    3.17e-15   3.4e-15  3.14e-15 ; ...
       128    9.46e-07    6.25e-15  7.04e-15  6.90e-15 ; ...
       256    7.22e-07    1.11e-14  1.37e-14  1.36e-14 ; ...
       512    2.11e-06    2.11e-14   2.8e-14  2.76e-14 ; ...
      1024    1.57e-06    4.12e-14  5.04e-14  5.20e-14 ; ...
      2048     4.2e-06    8.59e-14   1.1e-13  1.11e-13 ; ...
      4096    3.14e-06    1.65e-13  2.26e-13  2.25e-13 ; ...
      8192    9.02e-06    3.27e-13  4.29e-13  4.30e-13 ; ...
 1.638e+04    6.52e-06    6.66e-13  8.68e-13  8.68e-13 ; ...
 3.277e+04     1.9e-05     1.3e-12  1.69e-12  1.69e-12 ; ...
 6.554e+04    1.36e-05    2.63e-12  3.44e-12  3.44e-12 ; ...
 1.311e+05    3.85e-05    5.23e-12  6.94e-12  6.94e-12 ; ...
 2.621e+05    2.73e-05    1.03e-11  1.36e-11  1.36e-11 ; ...
 5.243e+05    7.15e-05    2.08e-11  2.78e-11  2.78e-11 ];

lam    = res(:,1);
SP_err = res(:,2);
DP_err = res(:,3);
CPU_err = res(:,4);
cmp_err = res(:,5);

figure(7)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[0.8 0.6]; set(gcf,'pos',pos);
loglog(lam,SP_err,'k-*',lam,DP_err,'k-.*',lam,CPU_err,'k-.o',lam,cmp_err,'k-.x')
xlabel('$\lambda$','interpreter','latex')
ylabel('$L_1$ error','interpreter','latex')
legend('poissinvf (GPU)','poissinv (GPU)','poissinv (CPU)','poisscinv (CPU)','Location','SouthEast')
axis([1 10^6 10^(-20) 10^(-3)])
print('-deps2','accuracy.eps')

end


%
% MATLAB function obtained from
% http://www.mathworks.co.uk/matlabcentral/fileexchange/
%  30347-sigplot/content/sigplot/sigplot/BasicFunctions/dsxy2figxy.m
%


function varargout = dsxy2figxy(varargin)
% dsxy2figxy -- Transform point or position from data space 
% coordinates into normalized figure coordinates 
% Transforms [x y] or [x y width height] vectors from data space
% coordinates to normalized figure coordinates in order to locate
% annotation objects within a figure. These objects are: arrow, 
% doublearrow, textarrow, ellipse line, rectangle, textbox 
%
% Syntax:
%    [figx figy] = dsxy2figxy([x1 y1],[x2 y2])  % GCA is used
%    figpos      = dsxy2figxy([x1 y1 width height])
%    [figx figy] = dsxy2figxy(axes_handle, [x1 y1],[x2 y2])
%    figpos      = dsxy2figxy(axes_handle, [x1 y1 width height])
%
% Usage: Obtain a position on a plot in data space and  
%        apply this function to locate an annotation there, e.g., 
%   [axx axy] = ginput(2); (input is in data space)
%   [figx figy] = dsxy2figxy(gca, axx, axy);  (now in figure space)
%   har = annotation('textarrow',figx,figy); 
%   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')']) 
%
%   Copyright 2006-2009 The MathWorks, Inc. 

% Obtain arguments (limited argument checking is done)
% Determine if axes handle is specified
if length(varargin{1}) == 1 && ishandle(varargin{1}) ...
                            && strcmp(get(varargin{1},'type'),'axes')	
	hAx = varargin{1};
	varargin = varargin(2:end); % Remove arg 1 (axes handle)
else
	hAx = gca;
end;

% Remaining args are either two point locations or a position vector
if length(varargin) == 1        % Assume a 4-element position vector
	pos = varargin{1};
else
	[x,y] = deal(varargin{:});  % Assume two pairs (start, end points)
end

% Get limits
axun = get(hAx,'Units');
set(hAx,'Units','normalized');  % Make axes units normalized 
axpos = get(hAx,'Position');    % Get axes position
axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
axwidth = diff(axlim(1:2));
axheight = diff(axlim(3:4));

% Transform from data space coordinates to normalized figure coordinates 
if exist('x','var')     % Transform a and return pair of points
	varargout{1} = (x - axlim(1)) * axpos(3) / axwidth + axpos(1);
	varargout{2} = (y - axlim(3)) * axpos(4) / axheight + axpos(2);
else                    % Transform and return a position rectangle
	pos(1) = (pos(1) - axlim(1)) / axwidth * axpos(3) + axpos(1);
	pos(2) = (pos(2) - axlim(3)) / axheight * axpos(4) + axpos(2);
	pos(3) = pos(3) * axpos(3) / axwidth;
	pos(4) = pos(4) * axpos(4 )/ axheight;
	varargout{1} = pos;
end

% Restore axes units
set(hAx,'Units',axun)

end
