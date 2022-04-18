classdef TB_patch_multi < TB_patch

properties
   
   MP    % MDTB-spline multi-patch
   H     % extraction matrix
   Hpol  % polynomial extraction matrix
   
end
   
methods

   function P = TB_patch_multi(MP, H, Hpol)
      
      % Construction of a TB-spline patch encapsulating a multi-patch
      %
      % INPUT
      %   MP    : MDTB-spline multi-patch
      %   H     : non-periodic extraction matrix
      %   Hpol  : non-periodic polynomial extraction matrix (optional)
      %
      % OUTPUT
      %   P     : TB-spline patch

      if nargin < 3
         Hpol = [];
      end
      n = size(H, 1);
      [a, b] = MDTB_domain(MP);
      P@TB_patch(n-1, [a, b]);
      P.p = [MP.P.p];
      P.MP = MP;
      P.H = H;
      P.Hpol = Hpol;
      
   end
   
   function gg = TB_greville(P)

      % Computation of Greville points (if possible)
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   gg    : vector of Greville points

      gg = MDTB_greville(P.MP, P.H);

   end
   
   function M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all MDTB-splines in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   xx    : vector of evaluation points
      %   cl    : closed domain if true (optional)
      %
      % OUTPUT
      %   M     : evaluation matrix 

      if nargin < 3
         M = MDTB_evaluation_all(P.MP, P.H, xx);
      else
         tol = 1e-12;
         M = zeros(P.n, length(xx));
         j = (xx >= P.U(1) & xx < P.U(end));
         if cl
            j(abs(xx-P.U(end))<tol) = true;
         end
         if any(j)
            M(:, j) = MDTB_evaluation_all(P.MP, P.H, xx(j));
         end
      end

   end
   
   function ss = TB_evaluation_spline(P, cc, xx)

      % Evaluation of a spline in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   cc    : vector of coefficients
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : vector of spline values

      ss = MDTB_evaluation_spline(P.MP, P.H, cc, xx);

   end  
   
   function ss = TB_evaluation_curve(P, cc, xx)

      % Evaluation of a spline curve in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   cc    : matrix of control points (1D, 2D or 3D)
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : matrix of spline curve points

      ss = MDTB_evaluation_curve(P.MP, P.H, cc, xx);

   end
   
   function K = TB_diffend_all(P, r, el)

      % Full differentiation of all MDTB-splines at one end point
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
      if el
         i = 1;
      else
         i = length(P.MP.P);
      end         
      K = P.H(:, P.MP.mu(i)+1:P.MP.mu(i+1)) * TB_diffend_all(P.MP.P(i), r, el);
      
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
         M = MDTB_differentiation_all(P.MP, P.H, r, xx);
      else
         tol = 1e-12;
         M = zeros(P.n, length(xx));
         j = (xx >= P.U(1) & xx < P.U(end));
         if cl
            j(abs(xx-P.U(end))<tol) = true;
         end
         if any(j)
            M(:, j) = MDTB_differentiation_all(P.MP, P.H, r, xx(j));
         end
      end
      
   end

   function ss = TB_differentiation_spline(P, r, cc, xx)

      % Differentiation of a spline in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   r     : order of derivative
      %   cc    : vector of coefficients
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : vector of r-th order derivative spline values

      ss = MDTB_differentiation_spline(P.MP, P.H, r, cc, xx);

   end
   
   function ss = TB_differentiation_curve(P, r, cc, xx)

      % Differentiation of a spline curve in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   r     : order of derivative
      %   cc    : matrix of control points (1D, 2D or 3D)
      %   xx    : vector of evaluation points
      %
      % OUTPUT
      %   ss    : matrix of r-th order derivative spline curve points

      ss = MDTB_differentiation_curve(P.MP, P.H, r, cc, xx);

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

      if isempty(P.Hpol)
         gg = TB_greville(P);
      else
         gg = MDTB_greville_poly(P.MP, P.Hpol);
      end

   end

end

end
