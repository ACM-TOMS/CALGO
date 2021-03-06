% Plots solution and grid levels from data files sol.dat and grid.dat
% generated by WRTUNI.f
% NB. pcolor with default shading colors a cell with the lowerleft value and
% ignores the last row and column.
%
nxb=input('nX base grid? ');
nyb=input('nY base grid? ');
nzb=input('nZ base grid? ');

load sol.dat
load grid.dat
[n,npde]=size(sol);

unilev=floor(log(n/(nxb*nyb*nzb))/(log(2)*3)+1);
lmul = 2^(unilev-1);
nx=nxb*2^(unilev-1); ny=nyb*2^(unilev-1); nz=nzb*2^(unilev-1);

while 1
   sldir = input('Enter slice coordinate direction, 1=x, 2=y, 3=z; 0 quits: ');
   if sldir == 0, break, end
   if sldir == 1
      slind = input(['Base grid index for slice, (0..',int2str(nxb),')? ']);
   elseif sldir == 2
      slind = input(['Base grid index for slice, (0..',int2str(nyb),')? ']);
   elseif sldir == 3
      slind = input(['Base grid index for slice, (0..',int2str(nzb),')? ']);
   end
   slind = slind*lmul;

   for ic=1:npde
      Umin=input(['min. sol. value comp. ',int2str(ic),'? ']);
      Umax=input(['max. sol. value comp. ',int2str(ic),'? ']);
      if sldir == 1
	 U = zeros(nz+1,ny+1);
	 for k = 0:nz
	 for j = 0:ny
            U(k+1,j+1) = sol(((ny+1)*k+j)*(nx+1)+slind+1,ic);
	 end
	 end
      elseif sldir == 2
	 U = zeros(nz+1,nx+1);
	 for k = 0:nz
	 for i = 0:nx
            U(k+1,i+1) = sol(((ny+1)*k+slind)*(nx+1)+i+1,ic);
	 end
	 end
      elseif sldir == 3
	 U = zeros(ny+1,nx+1);
	 for j = 0:ny
	 for i = 0:nx
            U(j+1,i+1) = sol(((ny+1)*slind+j)*(nx+1)+i+1,ic);
	 end
	 end
      end
      figure(ic);
      colormap(jet); pcolor(U); shading('interp');
      if Umin < Umax
         caxis([Umin Umax]);
      end
      keyboard
   end

   if sldir == 1
      G = zeros(nz+2,ny+2);
      for k = 0:nz
      for j = 0:ny
         G(k+1,j+1) = grid(((ny+1)*k+j)*(nx+1)+slind+1);
      end
      end
   elseif sldir == 2
      G = zeros(nz+2,nx+2);
      for k = 0:nz
      for i = 0:nx
         G(k+1,i+1) = grid(((ny+1)*k+slind)*(nx+1)+i+1);
      end
      end
   elseif sldir == 3
      G = zeros(ny+2,nx+2);
      for j = 0:ny
      for i = 0:nx
         G(j+1,i+1) = grid(((ny+1)*slind+j)*(nx+1)+i+1);
      end
      end
   end
   figure(npde+1);
   grc=[1 1 1; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 0; 1 0 1];
   colormap(grc(1:unilev,:));
   pcolor(G);
   caxis([0.99 unilev]);
   keyboard
end
