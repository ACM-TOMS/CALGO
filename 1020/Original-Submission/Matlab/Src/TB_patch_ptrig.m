classdef TB_patch_ptrig < TB_patch_poly
   
properties
   
   w     % TB-spline parameter
   ss    % vector of scaling coefficients
   
end

methods
   
   function P = TB_patch_ptrig(p, xx, w)
      
      % Construction of a polynomial trigonometric TB-spline patch of even degree p
      %
      % INPUT
      %   p     : TB-spline degree
      %   xx    : vector of end points
      %   w     : TB-spline parameter (optional)
      %
      % OUTPUT
      %   P     : TB-spline patch  
      
      if nargin < 3
         w = 1;
      end
      if mod(p, 2) 
         p = p + 1;
      end
      P@TB_patch_poly(p, xx);
      P.w = abs(w);
      cwl = 2 * cos(P.w * (P.U(end) - P.U(1)) / 2);
      P.ss = eye(P.n, 1);
      for q = 3:2:P.n
         i = q:-1:3;
         P.ss(i) = P.ss(i-2)  + cwl * P.ss(i-1) + P.ss(i);
         P.ss(2) = cwl + P.ss(2);
      end
      
   end
   
   function M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all polynomial trigonometric TB-splines in given points
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
      if P.w * (P.U(end) - P.U(1)) < tol
         M = TB_evaluation_all@TB_patch_poly(P, xx, cl);
      else
         M = zeros(P.n, length(xx));
         j = (xx >= P.U(1) & xx < P.U(end));
         if cl
            j(abs(xx-P.U(end))<tol) = true;
         end
         if any(j)
            sswl = (sin(P.w * (P.U(end) - P.U(1)) / 2))^2;
            cswl = 2 * cos(P.w * (P.U(end) - P.U(1)) / 2) / sswl;
            swx1 = sin(P.w * (xx(j) - P.U(1)) / 2);
            swx2 = sin(P.w * (P.U(end) - xx(j)) / 2);
            xl1 = swx1.^2 / sswl;
            xl2 = swx1 .* swx2 * cswl;
            xl3 = swx2.^2 / sswl;
            M(1, j) = 1;
            for q = 3:2:P.n
               i = q:-1:3;
               M(i, j) = bsxfun(@times, xl1, M(i-2, j)) ...
                       + bsxfun(@times, xl2, M(i-1, j)) ...
                       + bsxfun(@times, xl3, M(i, j));
               M(2, j) = xl2 .* M(1, j) + xl3 .* M(2, j);
               M(1, j) = xl3 .* M(1, j);
            end
         end
      end

   end

   function K = TB_diffend_all(P, r, el)

      % Full differentiation of all polynomial trigonometric TB-splines 
      % at one end point up to a given order
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
      tol = 1e-12;
      if P.w * (P.U(end) - P.U(1)) < tol
         K = TB_diffend_all@TB_patch_poly(P, r, el);
      else
         K = zeros(P.n, r+1);
         if el
            K(1, 1) = 1;
            if P.p > 0 && r > 0
               cwl = 2 * cos(P.w * (P.U(end) - P.U(1)) / 2);
               swl = 2 * sin(P.w * (P.U(end) - P.U(1)) / 2);
               k = min(r, P.p);
               i = (1:k)';
               is1 = i .* P.ss(i+1) ./ P.ss(i);
               is2 = (P.n - i) .* P.ss(i) ./ P.ss(i+1);
               is3 = (P.p/2 + 1 - [i; k+1]) * cwl;
               for q = 1:r
                  k = min(q, P.p);
                  K(1:k+1, q+1) = ([0; is1(1:k) .* K(1:k, q)] ...
                     - [is2(1:k) .* K(2:k+1, q); 0] ...
                     - is3(1:k+1) .* K(1:k+1, q)) * P.w / swl;
               end
            end
         else
            K(P.n, 1) = 1;
            if P.p > 0 && r > 0
               cwl = 2 * cos(P.w * (P.U(end) - P.U(1)) / 2);
               swl = 2 * sin(P.w * (P.U(end) - P.U(1)) / 2);
               k = min(r, P.p);            
               i = (P.p:-1:P.n-k)';
               is1 = i .* P.ss(i+1) ./ P.ss(i);
               is2 = (P.n - i) .* P.ss(i) ./ P.ss(i+1);
               is3 = (P.p/2 + 1 - [P.n; i]) * cwl;
               for q = 1:r
                  k = min(q, P.p);
                  K(P.n-k:P.n, q+1) = ([0; is1(k:-1:1) .* K(P.n-k:P.p, q)] ...
                     - [is2(k:-1:1) .* K(P.n-k+1:P.n, q); 0] ...
                     - is3(k+1:-1:1) .* K(P.n-k:P.n, q)) * P.w / swl;
               end
            end
         end
      end

   end

   function M = TB_differentiation_all(P, r, xx, cl)

      % Differentiation of all polynomial trigonometric TB-splines in given points
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
      tol = 1e-12;
      if P.w * (P.U(end) - P.U(1)) < tol
         M = TB_differentiation_all@TB_patch_poly(P, r, xx, cl);
      else
         M = TB_evaluation_all(P, xx, cl);
         if r > 0
            j = (xx >= P.U(1)-tol & xx <= P.U(end)+tol);
            if any(j)
               cwl = 2 * cos(P.w * (P.U(end) - P.U(1)) / 2);
               swl = 2 * sin(P.w * (P.U(end) - P.U(1)) / 2);
               i = (1:P.p)';
               is1 = i .* P.ss(i+1) ./ P.ss(i);
               is2 = (P.n - i) .* P.ss(i) ./ P.ss(i+1);
               is3 = (P.p/2 + 1 - [i; P.n]) * cwl;
               z = zeros(1, sum(j));
               for q = 1:r
                  M(:, j) = [z; bsxfun(@times, is1, M(1:end-1, j))] ...
                          - [bsxfun(@times, is2, M(2:end, j)); z] ...
                          - bsxfun(@times, is3, M(:, j));
               end
               M = M * (P.w / swl)^r;
            end
         end
      end

   end

end

end
