function errplotl(sol,eldata,ev,xy,x,y,fig)
%errplotl   plots solution and error estimate on L-shaped domain
%   errplotl(sol,eldata,ev,xy,x,y,fig)
%   input
%          sol        nodal solution vector 
%          eldata     element error vector
%          ev         element mapping matrix
%          xy         vertex coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          fig        figure number
%
%   IFISS function: DJS; 5 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
fprintf('plotting solution and estimated errors... ')
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
[II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
figure(fig)
subplot(221),contour(X,Y,xysol,20),axis('square')
axis('off'), ellx, title('Finite Element Solution')
subplot(222),mesh(X,Y,xysol),axis('square')
view(330,30)
%%
xx=xy(:,1); yy=xy(:,2);
nel=length(eldata);
% loop over elements    
for ielem = 1:nel
xl = xx(ev(ielem,:));
yl = yy(ev(ielem,:)); 
xc(ielem,1) = 0.25*sum(xl);
xc(ielem,2) = 0.25*sum(yl);
end
%
% interpolate to a cartesian product mesh
x=0.5*(x(1:end-1)+x(2:end));
y=0.5*(y(1:end-1)+y(2:end));
[X,Y]=meshgrid(x,y);
xysol = griddata(xc(:,1),xc(:,2),eldata,X,Y);
[II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
subplot(223),contour(X,Y,xysol,15),axis('square')
axis('off'), ellx, title('Estimated Error')
subplot(224),mesh(X,Y,xysol),axis('square')
view(330,30)
subplot(111)
fprintf('done\n')
return
