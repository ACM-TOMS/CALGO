function ddisp(w,beta,z,c)
%DDISP  Display results of Schwarz-Christoffel disk parameter problem.
%       DDISP(W,BETA,Z,C) displays the results of DPARAM in a pleasant
%       way.  
%
%	See also DPARAM, DPLOT.
%	
%	Written by Toby Driscoll.  Last updated 5/24/95.

disp(' ')
disp('          w               beta              z               arg(z)/pi')
disp(' -----------------------------------------------------------------------')
u = real(w);
v = imag(w);
x = real(z);
y = imag(z);
ang = angle(z)/pi;
ang(ang<=0) = ang(ang<=0) + 2;
for j = 1:length(w)
  if v(j) < 0
    s1 = '-';
  else
    s1 = '+';
  end
  if y(j) < 0
    s2 = '-';
  else
    s2 = '+';
  end
  disp(sprintf(' %8.5f %c %7.5fi    %8.5f   %8.5f %c %7.5fi    %14.12f',...
      u(j),s1,abs(v(j)),beta(j),x(j),s2,abs(y(j)),ang(j)));
  
end
disp(' ')
if imag(c) < 0
  s = '-';
else
  s = '+';
end
disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))
