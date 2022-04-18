% IM    interval maps:
%         maps xi in [-1,1] to x in [a,b] and
%         maps x in [a,b]  to xi in [-1,1]
% Called by:
%   1) grp.m, 2) frp.m, 3) inverseReprojection.m
% Last modified: October 17, 2007

function r = im(a,b,x,ch)
             
   switch ch
     case 0                          % maps xi in [-1,1] to x in [a,b]      x(xi)
      r = 0.5*(b-a).*x + 0.5*(b+a);
    case 1               
      r = (2.*x-b-a)./(b-a);         % maps x in [a,b]  to xi in [-1,1]    xi(x)
    otherwise
      r = x;
   end
