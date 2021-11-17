function display(p)
%   Show the polygon as a list of vertices in (x,y) form.
%   Copyright 1998 by Toby Driscoll.
%   $Id: display.m,v 2.2 2001/05/16 13:52:00 driscoll Exp $

w = p.vertex;
n = length(w);

fprintf('\n%s = \n\n',inputname(1))
if n==0
  fprintf('     []\n\n')
else
  x = real(w);
  y = imag(w);
  pm = '+-';
  s = 1 + (sign(y) < 0);
  for j = 1:n
    if ~isinf(w(j))
      if y(j)==0
	fprintf('  %7.4f\n',x(j))
      else
	fprintf('  %7.4f %c %6.4fi\n',x(j),pm(s(j)),abs(y(j)))
      end
    else
      fprintf('  %7s\n','Inf')
    end
  end
  fprintf('\n')
end
