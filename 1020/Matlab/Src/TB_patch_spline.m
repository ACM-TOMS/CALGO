classdef TB_patch_spline < TB_patch
   
methods

   function P = TB_patch_spline(p, xx, kk)
      
      % Construction of a polynomial B-spline patch of degree p with open knots
      %
      % INPUT
      %   p     : B-spline degree
      %   xx    : vector of break points
      %   kk    : smoothness vector (optional)
      %
      % OUTPUT
      %   P     : TB-spline patch
      %
      % If kk is a scalar, smoothness kk is imposed at break point xx(i+1),
      % if kk is a vector, smoothness kk(i) is imposed at break point xx(i+1),
      % for i = 1:length(xx)-2

      if nargin < 3
         kk = 0;
      end
      if length(kk) == 1
         kk = repmat(kk, 1, length(xx)-2);
      end
      U = repelem(xx, [p+1, p-kk, p+1]);
      P@TB_patch(p, U);
      P.n = length(U) - p - 1;
      
   end
   
   function gg = TB_greville(P)

      % Computation of polynomial spline Greville points
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   gg    : vector of Greville points

      gg = mean(P.U(hankel(2:P.p+1, P.p+1:P.p+P.n)), 1);

   end
   
   function M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all polynomial B-splines in given points
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
      for i = P.p+1:P.n
         j = (xx >= P.U(i) & xx < P.U(i+1));
         M(i, j) = 1;
         if cl && P.U(i+1) == P.U(end)
            M(i, abs(xx-P.U(end))<tol) = 1;
         end
      end
      for q = 1:P.p
         for i = P.p-q+1:P.n
            j = (xx >= P.U(i)-tol & xx <= P.U(i+q+1)+tol);
            l = P.U(i+q) - P.U(i);
            if l > 0 && any(j)
               M(i, j) = (xx(j) - P.U(i)) / l .* M(i, j);
            end
            l = P.U(i+q+1) - P.U(i+1);
            if l > 0 && any(j)
               M(i, j) = M(i, j) + (P.U(i+q+1) - xx(j)) / l .* M(i+1, j);
            end
         end
      end

   end
   
   function K = TB_diffend_all(P, r, el)

      % Full differentiation of all polynomial B-splines at one end point
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
            for q = P.p-r+1:P.p
               j = 1:r-P.p+q;
               for i = r-P.p+q+1:-1:1
                  l = P.U(i+P.p+1) - P.U(i+P.p+1-q);
                  if l > 0
                     K(i, j+1) = - q / l * K(i, j);
                  end
                  l = P.U(i+P.p) - P.U(i+P.p-q);
                  if i > 0 && l > 0
                     K(i, j+1) = K(i, j+1) + q / l * K(i-1, j);
                  end
               end
            end
         end
      else
         K(P.n, 1) = 1;
         if P.p > 0 && r > 0
            r = min(P.p, r);
            for q = P.p-r+1:P.p
               j = 1:r-P.p+q;
               for i = P.n-r+P.p-q:P.n
                  l = P.U(i+q) - P.U(i);
                  if l > 0
                     K(i, j+1) = q / l * K(i, j);
                  end
                  l = P.U(i+q+1) - P.U(i+1);
                  if i < P.n && l > 0
                     K(i, j+1) = K(i, j+1) - q / l * K(i+1, j);
                  end
               end
            end
         end
      end

   end

   function M = TB_differentiation_all(P, r, xx, cl)

      % Differentiation of all polynomial B-splines in given points
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
         Pr = P;
         Pr.p = P.p-r;
         Pr.n = P.n+r;
         Mr = TB_evaluation_all(Pr, xx, cl);
         for q = P.p-r+1:P.p
            for i = P.p-q+1:P.n
               j = (xx >= P.U(i)-tol & xx <= P.U(i+q+1)+tol);
               l = P.U(i+q) - P.U(i);
               if l > 0 && any(j)
                  Mr(i, j) = q / l * Mr(i, j);
               end
               l = P.U(i+q+1) - P.U(i+1);
               if l > 0 && any(j)
                  Mr(i, j) = Mr(i, j) - q / l * Mr(i+1, j);
               end
            end
         end
         M = Mr(1:P.n, :);
      else
         M = zeros(P.n, length(xx));
      end

   end

end

methods (Access = {?TB_patch, ?MDTB_patch})
   
   function gg = TB_greville_poly(P)

      % Computation of polynomial spline Greville points
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   gg    : vector of Greville points

      gg = TB_greville(P);

   end

end

end
