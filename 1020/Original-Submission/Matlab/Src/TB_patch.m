classdef (Abstract) TB_patch < matlab.mixin.Heterogeneous

properties
   
   p     % TB-spline degree
   n     % TB-spline dimension
   U     % vector of TB-spline points
   
end

methods
   
   function P = TB_patch(p, xx)
      
      % Construction of a TB-spline patch of degree p
      %
      % INPUT
      %   p     : TB-spline degree
      %   xx    : vector of end points
      %
      % OUTPUT
      %   P     : TB-spline patch  
      
      P.p = p;
      P.n = p+1;
      P.U = sort(xx);
      
   end
   
   function [a, b] = TB_domain(P)

      % Computation of end points of the domain
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   a     : left end point
      %   b     : right end point

      a = P.U(1);
      b = P.U(end);

   end
   
   function gg = TB_greville(P)

      % Computation of Greville points (if possible)
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   gg    : vector of Greville points

      xx = TB_greville_poly(P);
      M = TB_evaluation_all(P, xx);
      gg = xx / M;

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

      M = TB_evaluation_all(P, xx);
      ss = reshape(cc, 1, []) * M;

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

      M = TB_evaluation_all(P, xx);
      ss = cc * M;

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

      M = TB_differentiation_all(P, r, xx);
      ss = reshape(cc, 1, []) * M;

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

      M = TB_differentiation_all(P, r, xx);
      ss = cc * M;

   end

   function TB_visualization_all(P, n, varargin)

      % Visualization of all TB-splines
      %
      % INPUT
      %   P     : TB-spline patch
      %   n     : number of evaluation points (optional)
      %   specs : pass any number of plot specifications (optional)

      if nargin < 2
         n = 100;
      end
      [a, b] = TB_domain(P);
      xx = linspace(a, b, n);
      M = TB_evaluation_all(P, xx);
      plot(xx, M, varargin{:});

   end

   function TB_visualization_spline(P, cc, n, varargin)

      % Visualization of a spline
      %
      % INPUT
      %   P     : TB-spline patch
      %   cc    : vector of coefficients
      %   n     : number of evaluation points (optional)
      %   specs : pass any number of plot specifications (optional)

      if nargin < 3
         n = 100;
      end
      [a, b] = TB_domain(P);
      xx = linspace(a, b, n);
      ss = TB_evaluation_spline(P, cc, xx);
      plot(xx, ss, varargin{:});

   end

   function TB_visualization_curve(P, cc, n, varargin)

      % Visualization of a spline curve
      %
      % INPUT
      %   P     : TB-spline patch
      %   cc    : matrix of control points (1D, 2D or 3D)
      %   n     : number of evaluation points (optional)
      %   specs : pass any number of plot specifications (optional)

      if nargin < 3
         n = 100;
      end
      [a, b] = TB_domain(P);
      xx = linspace(a, b, n);
      ss = TB_evaluation_curve(P, cc, xx);
      switch size(cc, 1)
      case 1
         plot(xx, ss, varargin{:});
      case 2
         plot(ss(1, :), ss(2, :), varargin{:});
      case 3
         plot3(ss(1, :), ss(2, :), ss(3, :), varargin{:});
      end

   end

   function ccd = TB_conversion(Pd, Ps, ccs, sh)

      % Conversion from source to destination TB-spline form
      %
      % INPUT
      %   Pd    : destination TB-spline patch
      %   Ps    : source TB-spline patch
      %   ccs   : source coefficient vector
      %   sh    : shift of the source patch (optional)
      %
      % OUTPUT
      %   ccd   : destination coefficient vector

      if nargin < 4
         sh = 0;
      end
      gg = TB_greville_poly(Pd);
      Md = TB_evaluation_all(Pd, gg);
      sss = TB_evaluation_spline(Ps, ccs, gg - sh);
      ccd = sss / Md;

   end

end

methods (Access = {?TB_patch, ?MDTB_patch})
   
   function gg = TB_greville_poly(P)

      % Computation of polynomial Greville points
      %
      % INPUT
      %   P     : TB-spline patch
      %
      % OUTPUT
      %   gg    : vector of Greville points

      gg = (P.U(end) - P.U(1)) / P.p * (0:P.p) + P.U(1);

   end
   
end

methods (Abstract)
   
   M = TB_evaluation_all(P, xx, cl)

      % Evaluation of all TB-splines in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   xx    : vector of evaluation points
      %   cl    : closed domain if true (optional)
      %
      % OUTPUT
      %   M     : evaluation matrix 

   K = TB_diffend_all(P, r, el)

      % Full differentiation of all TB-splines at one end point
      % up to a given order
      %
      % INPUT
      %   P     : TB-spline patch
      %   r     : max order of derivative
      %   el    : left end if true, right end otherwise (optional)
      %
      % OUTPUT
      %   K     : differentiation matrix at end point up r-th order

   M = TB_differentiation_all(P, r, xx, cl)

      % Differentiation of all TB-splines in given points
      %
      % INPUT
      %   P     : TB-spline patch
      %   r     : order of derivative
      %   xx    : vector of evaluation points
      %   cl    : closed domain if true (optional)
      %
      % OUTPUT
      %   M     : differentiation matrix 

end

end
