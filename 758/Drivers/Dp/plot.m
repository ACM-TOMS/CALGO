% Plots solution and grid levels from data files sol.dat and grid.dat
% generated by WRTUNI.f
% NB. pcolor with default shading colors a cell with the lowerleft value and
% ignores the last row and column.
%
nxb=input('nX base grid? ');
nyb=input('nY base grid? ');

load sol.dat
load grid.dat
[n,npde]=size(sol);

unilev=floor(log(n/(nxb*nyb))/(log(2)*2)+1)
nx=nxb*2^(unilev-1); ny=nyb*2^(unilev-1);

for ic=1:npde
   Umin=input(['min. sol. value comp. ',int2str(ic),'? '])
   Umax=input(['max. sol. value comp. ',int2str(ic),'? '])
   U=zeros(ny+1,nx+1);
   for j = 0:ny
   for i = 0:nx
      U(j+1,i+1) = sol(j*(nx+1)+i+1,ic);
   end
   end
   figure(ic);
   colormap(jet); pcolor(U); shading('interp');
   if Umin < Umax
      caxis([Umin Umax])
   end;
   keyboard
end;

G=zeros(ny+2,nx+2);
for j = 0:ny
for i = 0:nx
   G(j+1,i+1) = grid(j*(nx+1)+i+1);
end
end
figure(npde+1);
grc=[1 1 1; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 0; 1 0 1];
colormap(grc(1:unilev,:));
pcolor(G);
caxis([0.99 unilev])
