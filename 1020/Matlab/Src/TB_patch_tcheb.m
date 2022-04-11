classdef TB_patch_tcheb < TB_patch

properties
   
   W     % TB-spline root parameters
   mu    % cummulative dimension
   C     % representation matrix
   
end
   
methods
   
   function P = TB_patch_tcheb(p, xx, ww, mm)
      
      % Construction of a TB-spline patch of degree p based on a linear
      % differential operator with constant coefficients
      %
      % INPUT
      %   p     : TB-spline degree
      %   xx    : vector of end points
      %   ww    : TB-spline roots (optional)
      %   mm    : TB-spline multiplicities (optional)
      %
      % OUTPUT
      %   P     : TB-spline patch  
      
      P@TB_patch(p, xx);
      if nargin < 3
         P.W = [0, 0, 0, p+1];
         P.mu = [0, p+1];
      else
         if nargin < 4
            mm = 1;
         end
         if length(mm) == 1
            mm = repmat(mm, 1, length(ww));
         end
         m = dot(logical(imag(ww)) + 1, mm);
         if P.n > m
            ww(end+1) = 0;
            mm(end+1) = P.n - m;
         elseif P.n < m
            P.n = m;
            P.p = m - 1;
         end
         i = mm >= 1;
         wa = reshape(ww(i), [], 1);
         [wu, ~, wi] = unique([real(wa), abs(imag(wa))], 'rows');    
         P.W = zeros(size(wu, 1), 4);
         P.W(:, 2:3) = wu;
         P.W(:, 4) = accumarray(wi, mm(i));
         P.W(logical(P.W(:, 2)), 1) = 1;
         P.W(logical(P.W(:, 3)), 1) = 2;
         P.W(all(P.W(:, 2:3), 2), 1) = 3;
         P.mu = [0, cumsum((logical(P.W(:, 3)) + 1) .* P.W(:, 4))'];
      end
      P.C = TB_representation_all(P);
      
   end
   
   function M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all linear differential TB-splines in given points
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
         xl = xx(j) - P.U(1);
         E = zeros(P.n, length(xl));
         m = max(P.W(:, 4)) - 1;
         X = [ones(size(xl)); cumprod((1 ./ (1:m))' * xl, 1)];
         for i = 1:size(P.W, 1)
            switch P.W(i, 1)
            case 0
               Ei = X(1:P.W(i, 4), :);
            case 1
               ewxl = exp(P.W(i, 2) * xl);
               Ei = bsxfun(@times, X(1:P.W(i, 4), :), ewxl);
            case 2
               cwxl = cos(P.W(i, 3) * xl);
               swxl = sin(P.W(i, 3) * xl);
               Ei = [bsxfun(@times, X(1:P.W(i, 4), :), cwxl);
                     bsxfun(@times, X(1:P.W(i, 4), :), swxl)];
            case 3
               ewxl = exp(P.W(i, 2) * xl);
               ecwxl = ewxl .* cos(P.W(i, 3) * xl);
               eswxl = ewxl .* sin(P.W(i, 3) * xl);
               Ei = [bsxfun(@times, X(1:P.W(i, 4), :), ecwxl);
                     bsxfun(@times, X(1:P.W(i, 4), :), eswxl)];
            end
            E(P.mu(i)+1:P.mu(i+1), :) = Ei;
         end
         M(:, j) = P.C * E;
      end
      
   end

   function K = TB_diffend_all(P, r, el)

      % Full differentiation of all linear differential TB-splines at one 
      % end point up to a given order
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
            E = zeros(P.n, r);
            for i = 1:size(P.W, 1)
               switch P.W(i, 1)
               case 0
                  Ei = [zeros(1, r); eye(P.W(i, 4)-1, r)];
               case 1
                  Ei = [repmat(P.W(i, 2), 1, r); eye(P.W(i, 4)-1, r)];
                  for k = 2:r
                     ki = min(k, P.W(i, 4)) - 1;
                     Ei(1:ki+1, k) = [0; Ei(1:ki, k-1)] + P.W(i, 2) * Ei(1:ki+1, k-1);
                  end
               case 2
                  Eci = [zeros(1, r); eye(P.W(i, 4)-1, r)];
                  Esi = [repmat(P.W(i, 3), 1, r); zeros(P.W(i, 4)-1, r)];
                  for k = 2:r
                     ki = min(k, P.W(i, 4)) - 1;
                     Eci(1:ki+1, k) = [0; Eci(1:ki, k-1)] - P.W(i, 3) * Esi(1:ki+1, k-1);
                     Esi(1:ki+1, k) = [0; Esi(1:ki, k-1)] + P.W(i, 3) * Eci(1:ki+1, k-1);
                  end
                  Ei = [Eci; Esi];
               case 3      
                  Eci = [repmat(P.W(i, 2), 1, r); eye(P.W(i, 4)-1, r)];
                  Esi = [repmat(P.W(i, 3), 1, r); zeros(P.W(i, 4)-1, r)];
                  for k = 2:r
                     ki = min(k, P.W(i, 4)) - 1;
                     Eci(1:ki+1, k) = [0; Eci(1:ki, k-1)] ...
                        + P.W(i, 2) * Eci(1:ki+1, k-1) - P.W(i, 3) * Esi(1:ki+1, k-1);
                     Esi(1:ki+1, k) = [0; Esi(1:ki, k-1)] ...
                        + P.W(i, 2) * Esi(1:ki+1, k-1) + P.W(i, 3) * Eci(1:ki+1, k-1);
                  end
                  Ei = [Eci; Esi];
               end
               E(P.mu(i)+1:P.mu(i+1), :) = Ei;
            end
            K(:, 2:end) = triu(P.C * E, -1);
         end
      else
         K(P.n, 1) = 1;
         if P.p > 0 && r > 0
            l = P.U(end) - P.U(1);
            E = zeros(P.n, r);
            m = max(P.W(:, 4)) - 1;
            X = [1; cumprod(l ./ (1:m)', 1)];
            for i = 1:size(P.W, 1)
               switch P.W(i, 1)
               case 0
                  if P.W(i, 4) == 1
                     Ei = zeros(1, r);
                  else
                     Ei = [zeros(1, r); toeplitz(X(1:P.W(i, 4)-1), eye(1, r))];
                  end
               case 1
                  ewl = exp(P.W(i, 2) * l);
                  Xewl = X(1:P.W(i, 4)) * ewl;     
                  Ei = zeros(P.W(i, 4), r);
                  Ei(:, 1) = [0; Xewl(1:end-1)] + P.W(i, 2) * Xewl;
                  for k = 2:r
                     Ei(:, k) = [0; Ei(1:end-1, k-1)] + P.W(i, 2) * Ei(:, k-1);
                  end
               case 2
                  cwl = cos(P.W(i, 3) * l);
                  swl = sin(P.W(i, 3) * l);
                  Xcwl = X(1:P.W(i, 4)) * cwl;     
                  Xswl = X(1:P.W(i, 4)) * swl;     
                  Eci = zeros(P.W(i, 4), r);
                  Esi = zeros(P.W(i, 4), r);
                  Eci(:, 1) = [0; Xcwl(1:end-1)] - P.W(i, 3) * Xswl;
                  Esi(:, 1) = [0; Xswl(1:end-1)] + P.W(i, 3) * Xcwl;
                  for k = 2:r
                     Eci(:, k) = [0; Eci(1:end-1, k-1)] - P.W(i, 3) * Esi(:, k-1);
                     Esi(:, k) = [0; Esi(1:end-1, k-1)] + P.W(i, 3) * Eci(:, k-1);
                  end
                  Ei = [Eci; Esi];
               case 3      
                  ecwl = exp(P.W(i, 2) * l) * cos(P.W(i, 3) * l);
                  eswl = exp(P.W(i, 2) * l) * sin(P.W(i, 3) * l);
                  Xcwl = X(1:P.W(i, 4)) * ecwl;     
                  Xswl = X(1:P.W(i, 4)) * eswl;     
                  Eci = zeros(P.W(i, 4), r);
                  Esi = zeros(P.W(i, 4), r);
                  Eci(:, 1) = [0; Xcwl(1:end-1)] + P.W(i, 2) * Xcwl - P.W(i, 3) * Xswl;
                  Esi(:, 1) = [0; Xswl(1:end-1)] + P.W(i, 2) * Xswl + P.W(i, 3) * Xcwl;
                  for k = 2:r
                     Eci(:, k) = [0; Eci(1:end-1, k-1)] ...
                        + P.W(i, 2) * Eci(:, k-1) - P.W(i, 3) * Esi(:, k-1);
                     Esi(:, k) = [0; Esi(1:end-1, k-1)] ...
                        + P.W(i, 2) * Esi(:, k-1) + P.W(i, 3) * Eci(:, k-1);
                  end
                  Ei = [Eci; Esi];
               end
               E(P.mu(i)+1:P.mu(i+1), :) = Ei;
            end
            i = P.n:-1:1;
            K(i, 2:end) = triu(P.C(i, :) * E, -1);
         end
      end

   end

   function M = TB_differentiation_all(P, r, xx, cl)

      % Differentiation of all linear differential TB-splines in given points
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
            xl = xx(j) - P.U(1);
            zxl = zeros(1, length(xl));
            E = zeros(P.n, length(xl));
            m = max(P.W(:, 4)) - 1;
            X = [ones(size(xl)); cumprod((1 ./ (1:m))' * xl, 1)];
            for i = 1:size(P.W, 1)
               switch P.W(i, 1)
               case 0
                  Ei = zeros(P.W(i, 4), length(xl));
                  Ei(r+1:end, :) = X(1:P.W(i, 4)-r, :);
               case 1
                  ewxl = exp(P.W(i, 2) * xl);
                  Ei = bsxfun(@times, X(1:P.W(i, 4), :), ewxl);
                  for k = 1:r
                     Ei = [zxl; Ei(1:end-1, :)] + P.W(i, 2) * Ei;
                  end
               case 2
                  cwxl = cos(P.W(i, 3) * xl);
                  swxl = sin(P.W(i, 3) * xl);
                  Eci = bsxfun(@times, X(1:P.W(i, 4), :), cwxl);
                  Esi = bsxfun(@times, X(1:P.W(i, 4), :), swxl);                     
                  for k = 1:r
                     Ezi = Eci;
                     Eci = [zxl; Eci(1:end-1, :)] - P.W(i, 3) * Esi;
                     Esi = [zxl; Esi(1:end-1, :)] + P.W(i, 3) * Ezi;
                  end
                  Ei = [Eci; Esi];
               case 3      
                  ewxl = exp(P.W(i, 2) * xl);
                  ecwxl = ewxl .* cos(P.W(i, 3) * xl);
                  eswxl = ewxl .* sin(P.W(i, 3) * xl);
                  Eci = bsxfun(@times, X(1:P.W(i, 4), :), ecwxl);
                  Esi = bsxfun(@times, X(1:P.W(i, 4), :), eswxl);
                  for k = 1:r
                     Ezi = Eci;
                     Eci = [zxl; Eci(1:end-1, :)] ...
                         + P.W(i, 2) * Eci - P.W(i, 3) * Esi;
                     Esi = [zxl; Esi(1:end-1, :)] ...
                         + P.W(i, 2) * Esi + P.W(i, 3) * Ezi;
                  end
                  Ei = [Eci; Esi];
               end
               E(P.mu(i)+1:P.mu(i+1), :) = Ei;
            end
            M(:, j) = P.C * E;
         end
      end

   end

end

methods (Access = protected)
   
   function C = TB_representation_all(P)

      % Representation of all linear differential TB-splines
      % (auxiliary function for TB_evaluation_all)
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   C     : representation matrix 

      l = P.U(end) - P.U(1);
      M0 = zeros(P.n);
      M1 = zeros(P.n);
      m = max(P.W(:, 4)) - 1;
      X = [1; cumprod(l ./ (1:m)', 1)];
      for i = 1:size(P.W, 1)
         switch P.W(i, 1)
         case 0
            M0i = eye(P.W(i, 4), P.n);
            M1i = toeplitz(X(1:P.W(i, 4)), eye(1, P.n));
         case 1
            ewl = exp(P.W(i, 2) * l);
            M0i = eye(P.W(i, 4), P.n);
            M1i = zeros(P.W(i, 4), P.n);
            M1i(:, 1) = X(1:P.W(i, 4)) * ewl;            
            for k = 1:P.p
               ki = min(k, P.W(i, 4)) - 1;
               M0i(1:ki+1, k+1) = [0; M0i(1:ki, k)] + P.W(i, 2) * M0i(1:ki+1, k);
               M1i(:, k+1) = [0; M1i(1:end-1, k)] + P.W(i, 2) * M1i(:, k);
            end
         case 2
            cwl = cos(P.W(i, 3) * l);
            swl = sin(P.W(i, 3) * l);
            M0ci = eye(P.W(i, 4), P.n);
            M0si = zeros(P.W(i, 4), P.n);
            M1ci = zeros(P.W(i, 4), P.n);
            M1si = zeros(P.W(i, 4), P.n);
            M1ci(:, 1) = X(1:P.W(i, 4)) * cwl;
            M1si(:, 1) = X(1:P.W(i, 4)) * swl;
            for k = 1:P.p
               ki = min(k, P.W(i, 4)) - 1;
               M0ci(1:ki+1, k+1) = [0; M0ci(1:ki, k)] - P.W(i, 3) * M0si(1:ki+1, k);
               M0si(1:ki+1, k+1) = [0; M0si(1:ki, k)] + P.W(i, 3) * M0ci(1:ki+1, k);
               M1ci(:, k+1) = [0; M1ci(1:end-1, k)] - P.W(i, 3) * M1si(:, k);
               M1si(:, k+1) = [0; M1si(1:end-1, k)] + P.W(i, 3) * M1ci(:, k);
            end
            M0i = [M0ci; M0si];
            M1i = [M1ci; M1si];
          case 3
            ecwl = exp(P.W(i, 2) * l) * cos(P.W(i, 3) * l);
            eswl = exp(P.W(i, 2) * l) * sin(P.W(i, 3) * l);
            M0ci = eye(P.W(i, 4), P.n);
            M0si = zeros(P.W(i, 4), P.n);
            M1ci = zeros(P.W(i, 4), P.n);
            M1si = zeros(P.W(i, 4), P.n);
            M1ci(:, 1) = X(1:P.W(i, 4)) * ecwl;
            M1si(:, 1) = X(1:P.W(i, 4)) * eswl;
            for k = 1:P.p
               ki = min(k, P.W(i, 4)) - 1;
               M0ci(1:ki+1, k+1) = [0; M0ci(1:ki, k)] ...
                  + P.W(i, 2) * M0ci(1:ki+1, k) - P.W(i, 3) * M0si(1:ki+1, k);
               M0si(1:ki+1, k+1) = [0; M0si(1:ki, k)] ...
                  + P.W(i, 2) * M0si(1:ki+1, k) + P.W(i, 3) * M0ci(1:ki+1, k);
               M1ci(:, k+1) = [0; M1ci(1:end-1, k)] ...
                  + P.W(i, 2) * M1ci(:, k) - P.W(i, 3) * M1si(:, k);
               M1si(:, k+1) = [0; M1si(1:end-1, k)] ...
                  + P.W(i, 2) * M1si(:, k) + P.W(i, 3) * M1ci(:, k);
            end
            M0i = [M0ci; M0si];
            M1i = [M1ci; M1si];
         end
         M0(P.mu(i)+1:P.mu(i+1), :) = M0i;
         M1(P.mu(i)+1:P.mu(i+1), :) = M1i;
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
