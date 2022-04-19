function MP = MDTB_patch_ppoly(pp, xx, tp, ww)

   % Construction of an MDTB-spline multi-patch from polynomial segments
   % (algebraic/exponential/trigonometric)
   %
   % INPUT
   %   pp    : vector of TB-spline degrees
   %   xx    : vector of break points
   %   tp    : vector of TB-spline types (optional)
   %           types: 0 = TB_patch_poly
   %                  1 = TB_patch_pexp
   %                  2 = TB_patch_ptrig
   %   ww    : vector of TB-spline parameters (optional)
   %
   % OUTPUT
   %   MP    : MDTB-spline multi-patch

   if nargin < 3
      tp = 0;
   end
   if nargin < 4
      ww = 1;
   end
   if length(pp) == 1
      pp = repmat(pp, 1, length(xx)-1);
   end
   if length(tp) == 1
      tp = repmat(tp, 1, length(xx)-1);
   end
   if length(ww) == 1
      ww = repmat(ww, 1, length(xx)-1);
   end
   m = length(pp);
   PP = repmat(TB_patch_poly(0, []), 1, m);
   for i = 1:m
      switch tp(i)
      case 0
         PP(i) = TB_patch_poly(pp(i), xx(i:i+1));
      case 1
         PP(i) = TB_patch_pexp(pp(i), xx(i:i+1), ww(i));
      case 2
         PP(i) = TB_patch_ptrig(pp(i), xx(i:i+1), ww(i));
      end
   end
   MP = MDTB_patch(PP);

end
