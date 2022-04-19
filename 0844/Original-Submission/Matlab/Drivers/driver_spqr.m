% This is a driver for testspqr.  Given m and n it invokes
% testspqr with combinatons of fullR, pivot, and cn.  The
% other values required by testspqr remain constant.  You can
% make different drivers by changing these values.

% Coded by G. W. (Pete) Stewart
% Jun 23 2004

m = 20; n = 10;
for t = 1:2
   if t==2
      temp = m;
      m = n;
      n = temp;
   end
   svalmin = 3;
   gappos = 6;
   gapfac = 1e-6;
   tol = 1e-4;
   maxcols = n;
   for fullR = [1,0]
      for pivot = [1,0]
         cn = 1;
         testspqr
         pause
         if fullR==0 & pivot==0
            cn = 0;
            testspqr
            pause
         end
      end
   end
end
