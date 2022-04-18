classdef TB_patch_poly < TB_patch
   
methods
   
   function P = TB_patch_poly(p, xx)
      
      % Construction of a polynomial TB-spline patch of degree p
      %
      % INPUT
      %   p     : polynomial degree
      %   xx    : vector of end points
      %
      % OUTPUT
      %   P     : TB-spline patch  
      
      P@TB_patch(p, xx);
      
   end
   
   function gg = TB_greville(P)

      % Computation of polynomial Greville points
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   gg    : vector of Greville points

      gg = TB_greville_poly(P);

   end
   
   function M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all polynomial TB-splines in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   xx    : vector of evaluation points
      %   cl    : closed domain if true (optional)
      %
      % OUTPUT
      %   M     : evaluation matrix 

      if nargin < 3
         cl = true;
      end
      tol = 1e-12;
      M = zeros(P.n, length(xx));
      j = (xx >= P.U(1) & xx < P.U(end));
      if cl
         j(abs(xx-P.U(end))<tol) = true;
      end
      if any(j)
         l = P.U(end) - P.U(1);
         xl1 = (xx(j) - P.U(1)) / l;
         xl2 = (P.U(end) - xx(j)) / l;
         M(1, j) = 1;
         for q = 2:P.n
            i = q:-1:2;
            M(i, j) = bsxfun(@times, xl1, M(i-1, j)) + bsxfun(@times, xl2, M(i, j));
            M(1, j) = xl2 .* M(1, j);
         end
      end

   end

   function K = TB_diffend_all(P, r, el)

      % Full differentiation of all polynomial TB-splines at one end point
      % up to a given order
      %
      % INPUT
      %   P     : TB-spline patch
      %   r     : max order of derivative
      %   el    : left end if true, right end otherwise (optional)
      %
      % OUTPUT
      %   K     : differentiation matrix at end point up r-th order

      if nargin < 3
         el = true;
      end
      K = zeros(P.n, r+1);
      if el
         K(1, 1) = 1;
         if P.p > 0 && r > 0
            r = min(P.p, r);
            l = P.U(end) - P.U(1);
            for q = P.n-r+1:P.n
               ql = (q-1) / l;
               j = 1:r-P.n+q;
               i = r-P.p+q:-1:2;
               K(i, j+1) = ql * (K(i-1, j) - K(i, j));
               K(1, j+1) = -ql * K(1, j);
            end
         end
      else
         K(P.n, 1) = 1;
         if P.p > 0 && r > 0
            r = min(P.p, r);
            l = P.U(end) - P.U(1);
            for q = P.n-r+1:P.n
               ql = (q-1) / l;
               j = 1:r-P.n+q;
               i = P.n-r+P.n-q:P.p;
               K(i, j+1) = ql * (K(i, j) - K(i+1, j));
               K(P.n, j+1) = ql * K(P.n, j);
            end
         end
      end

   end

   function M = TB_differentiation_all(P, r, xx, cl)

      % Differentiation of all polynomial TB-splines in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   r     : order of derivative
      %   xx    : vector of evaluation points
      %   cl    : closed domain if true (optional)
      %
      % OUTPUT
      %   M     : differentiation matrix 

      if nargin < 4
         cl = true;
      end
      if r <= 0
         M = TB_evaluation_all(P, xx, cl);
      elseif r <= P.p
         tol = 1e-12;
         M = zeros(P.n, length(xx));
         Pr = TB_patch_poly(P.p-r, P.U);
         M(1:P.n-r, :) = TB_evaluation_all(Pr, xx, cl);
         j = (xx >= P.U(1)-tol & xx <= P.U(end)+tol);
         if any(j)
            l = P.U(end) - P.U(1);
            for q = P.n-r+1:P.n
               ql = (q-1) / l;
               i = q:-1:2;
               M(i, j) = ql * (M(i-1, j) - M(i, j));
               M(1, j) = -ql * M(1, j);
            end
         end
      else
         M = zeros(P.n, length(xx));
      end

   end

end

end
