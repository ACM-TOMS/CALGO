classdef MDTB_patch

properties
   
   P     % vector of TB-spline patches
   mu    % cummulative dimension
   
end

methods
   
   function MP = MDTB_patch(PP)

      % Construction of an MDTB-spline multi-patch from B-spline segments 
      %
      % INPUT
      %   PP    : vector of TB-spline patches
      %
      % OUTPUT
      %   MP    : MDB-spline multi-patch

      MP.P  = PP;
      MP.mu = [0, cumsum([PP.n])];

   end

   function H = MDTB_extraction(MP, rr)

      % Computation of multi-degree spline extraction matrix
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   rr    : MDTB-spline smoothness vector (optional)
      %
      % OUTPUT
      %   H     : extraction matrix

      if nargin < 2
         rr = 0;
      end
      m = length(MP.P) - 1;
      if length(rr) == 1
         rr = repmat(rr, 1, m);
      end
      H = speye(MP.mu(end));
      for i = 1:m
         r = min([rr(i), MP.P(i).p(end), MP.P(i+1).p(1)]);
         K = sparse([TB_diffend_all(MP.P(i), r, false);
                    -TB_diffend_all(MP.P(i+1), r, true)]);
         L = H(:, MP.mu(i)+1:MP.mu(i+2)) * K;
         for j = 0:r
            Hbar = MDTB_patch.MDTB_nullspace(L(:, j+1));
            H = Hbar * H;
            L = Hbar * L;
         end
      end

   end

   function H = MDTB_extraction_periodic(MP, rr, rp)

      % Computation of multi-degree spline extraction matrix with periodicity
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   rr    : MDTB-spline smoothness vector (optional)
      %   rp    : periodicity smoothness (optional)
      %           the value should be less than half the dimension (floored) 
      %           of the related non-periodic MDB-spline space
      %
      % OUTPUT
      %   H     : extraction matrix

      if nargin < 2
         rr = 0;
      end
      H = MDTB_extraction(MP, rr);
      if nargin > 2
         m = length(MP.P);
         r = min([rp, MP.P(m).p(end), MP.P(1).p(1)]);
         if size(H, 1) >= 2*(r+1)
            Hper = circshift(H, r+1);   
            K = sparse([TB_diffend_all(MP.P(m), r, false);
                       -TB_diffend_all(MP.P(1), r, true)]);
            Lper = Hper(:, [MP.mu(m)+1:MP.mu(m+1) 1:MP.mu(2)]) * K;
            for j = 0:r
               Hbar = MDTB_patch.MDTB_nullspace(Lper(:, j+1));
               Hper = Hbar * Hper;
               Lper = Hbar * Lper;
            end
            H = Hper;
         end
      end

   end
   
   function Hl = MDTB_extraction_local(MP, H, ip)

      % Computation of the local extraction matrix on a patch
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   ip    : index of patch
      %
      % OUTPUT
      %   Hl    : local extraction matrix

      Hl = H(:, MP.mu(ip)+1:MP.mu(ip+1));

   end

   function [a, b] = MDTB_domain(MP)

      % Computation of end points of the domain
      %
      % INPUT
      %   MP    : MDB-spline multi-patch
      %
      % OUTPUT
      %   a     : left end point
      %   b     : right end point

      a = MP.P(1).U(1);
      b = a + sum(arrayfun(@(P) P.U(end) - P.U(1), MP.P));

   end
   
   function gg = MDTB_greville(MP, H)

      % Computation of multi-degree Greville points (if possible)
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %
      % OUTPUT
      %   gg    : vector of Greville points

      x = 0;
      m = length(MP.P) - 1;
      dd = zeros(1, MP.mu(end));
      for i = 1:m
         dd(MP.mu(i)+1:MP.mu(i+1)) = x + TB_greville(MP.P(i));
         x = x + MP.P(i).U(end) - MP.P(i+1).U(1);
      end
      dd(MP.mu(m+1)+1:MP.mu(end)) = x + TB_greville(MP.P(end));
      gg = dd / H;
      [a, b] = MDTB_domain(MP);
      gg(gg < a) = a;
      gg(gg > b) = b;

   end

   function M = MDTB_evaluation_all(MP, H, xx)

      % Evaluation of all MDTB-splines in given points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   M     : evaluation matrix 

      yy = xx;
      tol = 1e-12;
      m = length(MP.P) - 1;
      N = zeros(MP.mu(end), length(xx));
      for i = 1:m
         j = (yy >= MP.P(i).U(1) & yy < MP.P(i).U(end));
         N(MP.mu(i)+1:MP.mu(i+1), j) = TB_evaluation_all(MP.P(i), yy(j), false);
         yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
      end
      j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end)+tol);
      N(MP.mu(m+1)+1:MP.mu(end), j) = TB_evaluation_all(MP.P(end), yy(j), true);
      M = H * N;
   
   end

   function ss = MDTB_evaluation_spline(MP, H, cc, xx)

      % Evaluation of a multi-degree spline in given points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   cc    : vector of coefficients
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : vector of spline values

      yy = xx;
      tol = 1e-12;
      m = length(MP.P) - 1;
      dd = reshape(cc, 1, []) * H;
      ss = zeros(1, length(xx));
      for i = 1:m
         dd_loc = dd(MP.mu(i)+1:MP.mu(i+1));
         j = (yy >= MP.P(i).U(1) & yy < MP.P(i).U(end));
         ss(j) = TB_evaluation_spline(MP.P(i), dd_loc, yy(j));
         yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
      end
      dd_loc = dd(MP.mu(m+1)+1:MP.mu(end));
      j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end)+tol);
      ss(j) = TB_evaluation_spline(MP.P(end), dd_loc, yy(j));

   end

   function ss = MDTB_evaluation_curve(MP, H, cc, xx)

      % Evaluation of a multi-degree spline curve in given points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   cc    : matrix of control points (1D, 2D or 3D)
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : matrix of spline curve points

      yy = xx;
      tol = 1e-12;
      m = length(MP.P) - 1;
      dd = cc * H;
      ss = zeros(size(cc, 1), length(xx));
      for i = 1:m
         dd_loc = dd(:, MP.mu(i)+1:MP.mu(i+1));
         j = (yy >= MP.P(i).U(1) & yy < MP.P(i).U(end));
         ss(:, j) = TB_evaluation_curve(MP.P(i), dd_loc, yy(j));
         yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
      end
      dd_loc = dd(:, MP.mu(m+1)+1:MP.mu(end));
      j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end)+tol);
      ss(:, j) = TB_evaluation_curve(MP.P(end), dd_loc, yy(j));

   end

   function M = MDTB_differentiation_all(MP, H, r, xx)

      % Differentiation of all MDTB-splines in given points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   r     : order of derivative
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   M     : differentiation matrix 

      yy = xx;
      tol = 1e-12;
      m = length(MP.P) - 1;
      N = zeros(MP.mu(end), length(xx));
      for i = 1:m
         j = (yy >= MP.P(i).U(1) & yy < MP.P(i).U(end));
         N(MP.mu(i)+1:MP.mu(i+1), j) = TB_differentiation_all(MP.P(i), r, yy(j), false);
         yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
      end
      j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end)+tol);
      N(MP.mu(m+1)+1:MP.mu(end), j) = TB_differentiation_all(MP.P(end), r, yy(j), true);
      M = H * N;

   end

   function ss = MDTB_differentiation_spline(MP, H, r, cc, xx)

      % Differentiation of a multi-degree spline in given points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   r     : order of derivative
      %   cc    : vector of coefficients
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : vector of r-th derivative spline values

      yy = xx;
      tol = 1e-12;
      m = length(MP.P) - 1;
      dd = reshape(cc, 1, []) * H;
      ss = zeros(1, length(xx));
      for i = 1:m
         dd_loc = dd(MP.mu(i)+1:MP.mu(i+1));
         j = (yy >= MP.P(i).U(1) & yy < MP.P(i).U(end));
         ss(j) = TB_differentiation_spline(MP.P(i), r, dd_loc, yy(j));
         yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
      end
      dd_loc = dd(MP.mu(m+1)+1:MP.mu(end));
      j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end)+tol);
      ss(j) = TB_differentiation_spline(MP.P(end), r, dd_loc, yy(j));

   end

   function ss = MDTB_differentiation_curve(MP, H, r, cc, xx)

      % Differentiation of a multi-degree spline curve in given points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   r     : order of derivative
      %   cc    : matrix of control points (1D, 2D or 3D)
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : matrix of r-th derivative spline curve points

      yy = xx;
      tol = 1e-12;
      m = length(MP.P) - 1;
      dd = cc * H;
      ss = zeros(size(cc, 1), length(xx));
      for i = 1:m
         dd_loc = dd(:, MP.mu(i)+1:MP.mu(i+1));
         j = (yy >= MP.P(i).U(1) & yy < MP.P(i).U(end));
         ss(:, j) = TB_differentiation_curve(MP.P(i), r, dd_loc, yy(j));
         yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
      end
      dd_loc = dd(:, MP.mu(m+1)+1:MP.mu(end));
      j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end)+tol);
      ss(:, j) = TB_differentiation_curve(MP.P(end), r, dd_loc, yy(j));

   end

   function MDTB_visualization_all(MP, H, n, varargin)

      % Visualization of all MDTB-splines
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   n     : number of evaluation points (optional)
      %   specs : pass any number of plot specifications (optional)

      if nargin < 3
         n = 100;
      end
      [a, b] = MDTB_domain(MP);
      xx = linspace(a, b, n);
      M = MDTB_evaluation_all(MP, H, xx);
      plot(xx, M, varargin{:});

   end
   
   function MDTB_visualization_spline(MP, H, cc, n, varargin)

      % Visualization of a multi-degree spline
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   cc    : vector of coefficients
      %   n     : number of evaluation points (optional)
      %   specs : pass any number of plot specifications (optional)

      if nargin < 4
         n = 100;
      end
      [a, b] = MDTB_domain(MP);
      xx = linspace(a, b, n);
      ss = MDTB_evaluation_spline(MP, H, cc, xx);
      plot(xx, ss, varargin{:});

   end
   
   function MDTB_visualization_curve(MP, H, cc, n, varargin)

      % Visualization of a multi-degree spline curve
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : extraction matrix
      %   cc    : matrix of control points (1D, 2D or 3D)
      %   n     : number of evaluation points (optional)
      %   specs : pass any number of plot specifications (optional)

      if nargin < 4
         n = 100;
      end
      [a, b] = MDTB_domain(MP);
      xx = linspace(a, b, n);
      ss = MDTB_evaluation_curve(MP, H, cc, xx);
      switch size(cc, 1)
      case 1
         plot(xx, ss, varargin{:});
      case 2
         plot(ss(1, :), ss(2, :), varargin{:});
      case 3
         plot3(ss(1, :), ss(2, :), ss(3, :), varargin{:});
      end

   end

   function ccd = MDTB_conversion(MPd, Hd, MPs, Hs, ccs, sh, gl)

      % Conversion from source to destination MDTB-spline form
      %
      % INPUT
      %   MPd   : destination MDTB-spline multi-patch
      %   Hd    : destination extraction matrix
      %   MPs   : source MDTB-spline multi-patch
      %   Hs    : source extraction matrix
      %   ccs   : source coefficient vector
      %   sh    : shift of the source patch (optional)
      %   gl    : global conversion if true (optional)
      %
      % OUTPUT
      %   ccd   : destination coefficient vector

      if nargin < 6
         sh = 0;
      end
      if nargin < 7
         gl = true;
      end
      if gl
         gg = MDTB_greville(MPd, Hd);
         Md = MDTB_evaluation_all(MPd, Hd, gg);
         sss = MDTB_evaluation_spline(MPs, Hs, ccs, gg - sh);
         ccd = sss / Md;
      else
         m = length(MPd.P);
         dds = reshape(ccs, 1, []) * Hs;
         ddd = zeros(1, MPd.mu(end));
         sh = sh + MPs.P(1).U(1) - MPd.P(1).U(1);
         for i = 1:m
            sh = sh - MPs.P(i).U(1) + MPd.P(i).U(1);
            dds_loc = dds(MPs.mu(i)+1:MPs.mu(i+1));
            ddd_loc = TB_conversion(MPd.P(i), MPs.P(i), dds_loc, sh);
            ddd(MPd.mu(i)+1:MPd.mu(i+1)) = ddd_loc;
            sh = sh + MPs.P(i).U(end) - MPd.P(i).U(end);
         end
         ccd = ddd / Hd;
      end

   end

end

methods (Access = {?TB_patch, ?MDTB_patch})
   
   function gg = MDTB_greville_poly(MP, Hpol)

      % Computation of polynomial multi-degree Greville points
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   Hpol  : polynomial extraction matrix
      %
      % OUTPUT
      %   gg    : vector of Greville points

      x = 0;
      m = length(MP.P) - 1;
      dd = zeros(1, MP.mu(end));
      for i = 1:m
         dd(MP.mu(i)+1:MP.mu(i+1)) = x + TB_greville_poly(MP.P(i));
         x = x + MP.P(i).U(end) - MP.P(i+1).U(1);
      end
      dd(MP.mu(m+1)+1:MP.mu(end)) = x + TB_greville_poly(MP.P(end));
      gg = dd / Hpol;
      [a, b] = MDTB_domain(MP);
      gg(gg < a) = a;
      gg(gg > b) = b;

   end

end

methods (Static, Access = protected)

   function Hbar = MDTB_nullspace(ll)

      % Computation of left null-space of a column of matrix L
      % (auxiliary function for MDB_extraction)
      %
      % INPUT
      %   ll    : a column of L
      %   
      % OUTPUT
      %   Hbar  : null-space matrix of ll

      q = length(ll); 
      i1 = find(ll, 1, 'first'); 
      i2 = find(ll, 1, 'last');
      dd = zeros(q-1, 2);
      dd(1:i1, 1) = 1;
      for j = i1:i2-2
         dd(j, 2) = -ll(j) / ll(j+1) * dd(j, 1);
         dd(j+1, 1) = 1 - dd(j, 2);
      end
      dd(i2-1:q-1, 2) = 1;
      Hbar = spdiags(dd, [0 1], q-1, q);

   end

end

end
