function solsurf(sol,xy,x,y,fig)
%solsurf   plots solution surface on square-shaped domain
%   solsurf(sol,xy,x,y,fig);
%   input
%          sol        nodal solution vector
%          xy         nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          fig        figure number
%
%   IFISS function: DJS; 4 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
fprintf('surfing solution... ')
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
figure(fig)
subplot(122),surf(X,Y,xysol),axis('square')
view(330,30)
shading interp
fprintf('done\n')
return