classdef TB_patch_gtrig < TB_patch

properties
   
   w     % TB-spline parameter
   t     % representation type
   m     % representation parameter
   C     % representation matrix
   
end

methods
   
   function P = TB_patch_gtrig(p, xx, w, t, m)
      
      % Construction of a generalized trigonometric TB-spline patch of degree p
      %
      % INPUT
      %   p     : TB-spline degree
      %   xx    : vector of end points
      %   w     : TB-spline parameter (optional)
      %   t     : representation type (optional)
      %   m     : representation parameter (optional)
      %
      % OUTPUT
      %   P     : TB-spline patch  
      
      if nargin < 3
         w = 1;
      end
      if nargin < 4
         t = abs(w) * (xx(end) - xx(1)) >= 3;
      end
      if nargin < 5
         m = 10;
      end
      P@TB_patch(p, xx);
      P.w = abs(w);
      P.t = t;
      P.m = m;
      P.C = TB_representation_all(P);
      
   end
   
   function M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all generalized trigonometric TB-splines in given points
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
         if P.t
            wxl = P.w * (xx(j) - P.U(1));
            E = [ones(size(wxl)); cumprod((1 ./ (1:P.p))' * wxl, 1)];
            E(P.p, :) = cos(wxl);   
            E(P.n, :) = sin(wxl);
            M(:, j) = P.C * E;
         else
            xl = xx(j) - P.U(1);
            ww = [1, cumprod(repmat(-P.w * P.w, 1, P.m), 2)];
            Ef = [ones(size(xl)); cumprod((1 ./ (1:P.p+2*P.m))' * xl, 1)];
            E = Ef(1:P.n, :);
            E(P.p, :) = ww * Ef(P.p:2:end, :);
            E(P.n, :) = ww * Ef(P.n:2:end, :);
            M(:, j) = P.C * E;
         end
      end

   end

   function K = TB_diffend_all(P, r, el)

      % Full differentiation of all generalized trigonometric TB-splines 
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
      K = zeros(P.n, r+1);
      if el
         K(1, 1) = 1;
         if P.t
            ww = cumprod(repmat(P.w, 1, r), 2);
            if P.p == 1 && r > 0
               E = zeros(2, r);
               E(1, 2:4:end) = -ww(2:4:end);
               E(1, 4:4:end) = ww(4:4:end);
               E(2, 1:4:end) = ww(1:4:end);
               E(2, 3:4:end) = -ww(3:4:end);
               K(:, 2:end) = P.C * E;
            elseif P.p > 1 && r > 0
               k = min(r, P.p-2);
               E = zeros(P.p, r);
               E(1:k, 1:k) = diag(ww(1:k));
               E(P.p-1, 2:4:end) = -ww(2:4:end);
               E(P.p-1, 4:4:end) = ww(4:4:end);
               E(P.p, 1:4:end) = ww(1:4:end);
               E(P.p, 3:4:end) = -ww(3:4:end);
               K(:, 2:end) = triu(P.C(:, 2:end) * E, -1);
            end
         else
            if P.p == 1 && r > 0
               E = zeros(2, r);
               E(2, 1) = 1;
               if r > 1
                  ww = cumprod(repmat(-P.w * P.w, 1, P.m), 2);
                  k = min(r, 1+2*P.m);
                  E(1, 2:2:k) = ww(1:ceil((k-1)/2));
                  k = min(r, 2+2*P.m);
                  E(2, 3:2:k) = ww(1:ceil((k-2)/2));
               end
               K(:, 2:end) = P.C * E;
            elseif P.p > 1 && r > 0
               E = eye(P.p, r);
               if r > P.p
                  ww = cumprod(repmat(-P.w * P.w, 1, P.m), 2);
                  k = min(r, P.p+2*P.m);
                  E(P.p-1, P.p+1:2:k) = ww(1:ceil((k-P.p)/2));
                  k = min(r, P.p+1+2*P.m);
                  E(P.p, P.p+2:2:k) = ww(1:ceil((k-P.p-1)/2));
               end
               K(:, 2:end) = triu(P.C(:, 2:end) * E, -1);
            end
         end
      else
         K(P.n, 1) = 1;
         if P.t
            ww = cumprod(repmat(P.w, 1, r), 2);
            wl = P.w * (P.U(end) - P.U(1));
            cwl = cos(wl);
            swl = sin(wl);
            if P.p == 1 && r > 0
               E = zeros(2, r);
               E(1:2, 1:4:end) = [-swl; cwl] * ww(1:4:end);
               E(1:2, 2:4:end) = [-cwl; -swl] * ww(2:4:end);
               E(1:2, 3:4:end) = [swl; -cwl] * ww(3:4:end);
               E(1:2, 4:4:end) = [cwl; swl] * ww(4:4:end);
               K(:, 2:end) = P.C * E;
            elseif P.p > 1 && r > 0
               k = min(r, P.p-2);
               E = zeros(P.p, r);
               if k > 0
                  E(:, 1:k) = toeplitz([1, cumprod(wl ./ (1:P.p-1), 2)], eye(1, k));
               end
               E(P.p-1, 1:4:end) = -swl;
               E(P.p-1, 2:4:end) = -cwl;
               E(P.p-1, 3:4:end) = swl;
               E(P.p-1, 4:4:end) = cwl;
               E(P.p, 1:4:end) = cwl;
               E(P.p, 2:4:end) = -swl;
               E(P.p, 3:4:end) = -cwl;
               E(P.p, 4:4:end) = swl;
               i = P.n:-1:1;
               K(i, 2:end) = triu(P.C(i, 2:end) * bsxfun(@times, E, ww), -1);
            end
         else
            l = P.U(end) - P.U(1);
            if P.p == 1 && r > 0
               ww = [1, cumprod(repmat(-P.w * P.w, 1, P.m), 2)];
               Ef = toeplitz([1, cumprod(l ./ (1:2*P.m), 2)], eye(1, r));
               E = [ww(2:end) * Ef(2:2:end, :); ww * Ef(1:2:end, :)];
               K(:, 2:end) = P.C * E;
            elseif P.p > 1 && r > 0
               ww = [1, cumprod(repmat(-P.w * P.w, 1, P.m), 2)];
               Ef = toeplitz([1, cumprod(l ./ (1:P.p-1+2*P.m), 2)], eye(1, r));
               E = Ef(1:P.p, :);
               E(P.p-1, :) = ww * Ef(P.p-1:2:end, :);
               E(P.p, :) = ww * Ef(P.p:2:end, :);
               i = P.n:-1:1;
               K(i, 2:end) = triu(P.C(i, 2:end) * E, -1);
            end
         end
      end

   end

   function M = TB_differentiation_all(P, r, xx, cl)

      % Differentiation of all generalized trigonometric TB-splines in given points
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
      else
         tol = 1e-12;
         M = zeros(P.n, length(xx));
         j = (xx >= P.U(1) & xx < P.U(end));
         if cl
            j(abs(xx-P.U(end))<tol) = true;
         end
         if any(j)
            k = min(r, P.p-1);
            if P.t
               wxl = P.w * (xx(j) - P.U(1));
               E = cumprod([ones(size(wxl)); (1 ./ (1:P.p-k))' * wxl], 1);
               switch mod(r, 4)
               case 0
                  E(P.p-k, :) = cos(wxl);
                  E(P.n-k, :) = sin(wxl);
               case 1
                  E(P.p-k, :) = -sin(wxl);
                  E(P.n-k, :) = cos(wxl);
               case 2
                  E(P.p-k, :) = -cos(wxl);
                  E(P.n-k, :) = -sin(wxl);
               case 3
                  E(P.p-k, :) = sin(wxl);
                  E(P.n-k, :) = -cos(wxl);
               end
               M(:, j) = P.w^r * (P.C(:, k+1:end) * E);
            else
               xl = xx(j) - P.U(1);
               ww = [1, cumprod(repmat(-P.w * P.w, 1, P.m), 2)];
               Ef = [ones(size(xl)); cumprod((1 ./ (1:P.p-k+2*P.m))' * xl, 1)];
               E = Ef(1:P.n-k, :);
               E(P.p-k, :) = ww * Ef(P.p-k:2:end, :);
               E(P.n-k, :) = ww * Ef(P.n-k:2:end, :);
               M(:, j) = P.C(:, k+1:end) * E;
            end
         end
      end

   end

end

methods (Access = protected)
   
   function C = TB_representation_all(P)

      % Representation of all generalized trigonometric TB-splines
      % (auxiliary function for TB_evaluation_all)
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   C     : representation matrix 

      if P.t
         wl = P.w * (P.U(end) - P.U(1));
         cwl = cos(wl);
         swl = sin(wl);
         M0 = eye(P.n);
         M0([P.p, P.n], [P.p, P.n]) = 0;
         M0(P.p, 1:4:end) = 1;
         M0(P.p, 3:4:end) = -1;
         M0(P.n, 2:4:end) = 1;
         M0(P.n, 4:4:end) = -1;
         M1 = toeplitz([1, cumprod(wl ./ (1:P.p), 2)], eye(1, P.n));
         M1(P.p, 1:4:end) = cwl;
         M1(P.p, 2:4:end) = -swl;
         M1(P.p, 3:4:end) = -cwl;
         M1(P.p, 4:4:end) = swl;
         M1(P.n, 1:4:end) = swl;
         M1(P.n, 2:4:end) = cwl;
         M1(P.n, 3:4:end) = -swl;
         M1(P.n, 4:4:end) = -cwl;
      else
         M0 = eye(P.n);
         l = P.U(end) - P.U(1);
         ww = [1, cumprod(repmat(-P.w * P.w, 1, P.m), 2)];
         M = toeplitz([1, cumprod(l ./ (1:P.p+2*P.m), 2)], eye(1, P.n));
         M1 = M(1:P.n, :);
         M1(P.p, :) = ww * M(P.p:2:end, :);
         M1(P.n, :) = ww * M(P.n:2:end, :);
      end
      cs = zeros(1, P.n);
      C = zeros(P.n);
      C(P.n, :) =  eye(1, P.n) / [M1(:, 1), M0(:, 1:P.p)];
      for i = 2:P.p
         cs = cs + C(P.n-i+2, :);
         cc = zeros(1, P.n);
         cc(i) = -cs * M1(:, i);
         C(P.n-i+1, :) = cc / [M1(:, 1:i), M0(:, 1:P.n-i)];
      end
      C(1, :) = eye(1, P.n) / [M0(:, 1), M1(:, 1:P.p)];

   end

end

end
