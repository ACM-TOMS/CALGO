function MP = MDTB_patch_tcheb(pp, xx, ww, mm)

   % Construction of an MDTB-spline multi-patch from Tchebycheff segments
   % based on linear differential operators with constant coefficients
   %
   % INPUT
   %   pp    : vector of TB-spline degrees
   %   xx    : vector of break points
   %   ww    : cell array of TB-spline roots (optional)
   %   mm    : cell array of TB-spline multiplicities (optional)
   %
   % OUTPUT
   %   MP    : MDTB-spline multi-patch

   if nargin < 3
      ww = {0};
   end
   if nargin < 4
      mm = {1};
   end
   if length(pp) == 1
      pp = repmat(pp, 1, length(xx)-1);
   end
   if length(ww) == 1
      ww = repmat(ww, 1, length(xx)-1);
   end
   if length(mm) == 1
      mm = repmat(mm, 1, length(xx)-1);
   end
   m = length(pp);
   PP = repmat(TB_patch_poly(0, []), 1, m);
   for i = 1:m
      PP(i) = TB_patch_tcheb(pp(i), xx(i:i+1), ww{i}, mm{i});
   end
   MP = MDTB_patch(PP);

end
