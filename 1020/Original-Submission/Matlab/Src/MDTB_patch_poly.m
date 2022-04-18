function [MP, rr] = MDTB_patch_poly(pp, xx, kk, mg)

   % Construction of an MDTB-spline multi-patch from polynomial segments
   % where consecutive polynomial segments of the same degree can be merged 
   % into a single TB-spline patch
   %
   % INPUT
   %   pp    : vector of polynomial degrees
   %   xx    : vector of break points
   %   kk    : smoothness vector (optional)
   %   mg    : same degree merged if true (optional)
   %
   % OUTPUT
   %   MP    : MDTB-spline multi-patch
   %   rr    : MDTB-spline smoothness vector (optional)
   %
   % If kk is a scalar, smoothness kk is imposed at break point xx(i+1),
   % if kk is a vector, smoothness kk(i) is imposed at break point xx(i+1),
   % for i = 1:length(xx)-2

   if nargin < 3
      kk = 0;
   end
   if nargin < 4
      mg = true;
   end
   if length(pp) == 1
      pp = repmat(pp, 1, length(xx)-1);
   end
   if length(kk) == 1
      kk = repmat(kk, 1, length(xx)-2);
   end
   if mg
      jj = [1 find(diff(pp))+1 length(pp)+1];
      m = length(jj) - 1;
      PP = repmat(TB_patch_poly(0, []), 1, m);
      for i = 1:m
         PP(i) = TB_patch_spline(pp(jj(i)), xx(jj(i):jj(i+1)), kk(jj(i):jj(i+1)-2));
      end
      MP = MDTB_patch(PP);
      rr = kk(jj(2:m)-1);
   else
      m = length(pp);
      PP = repmat(TB_patch_poly(0, []), 1, m);
      for i = 1:m
         PP(i) = TB_patch_poly(pp(i), xx(i:i+1));
      end
      MP = MDTB_patch(PP);
      rr = kk;
   end

end
