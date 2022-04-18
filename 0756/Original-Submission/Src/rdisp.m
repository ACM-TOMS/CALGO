function rdisp(w,beta,z,c,L)
%RDISP Display results of Schwarz-Christoffel rectangle parameter problem.
%       RDISP(W,BETA,RECT,Z,C) displays the results of RPARAM in a 
%       pleasant way.  
%
%	See also RPARAM, RPLOT.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
% Deduce corner locations
left = abs(real(z)-min(real(z))) < eps;
right = abs(real(z)-max(real(z))) < eps;
top = abs(imag(z)-max(imag(z))) < eps;
bot = abs(imag(z)-min(imag(z))) < eps;
corners = find(left+right+top+bot - 1);
c1 = find(abs(z-max(real(z))) < eps);
offset = find(corners==c1);
corners = corners([offset:4,1:offset-1]);
rect = z(corners);

disp(' ')
disp(' cnr          w               beta                      z              ')
disp(' ------------------------------------------------------------------------')
u = real(w);
v = imag(w);
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  cnr = find(j==corners);
  if isempty(cnr)
    cstr = '    ';
  else
    cstr = sprintf('  %i ',cnr);
  end
  if ~imag(z(j))
    disp(sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e',...
      cstr,u(j),s,abs(v(j)),beta(j),z(j)));
  else
    disp(sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e + %14.8ei',...
      cstr,u(j),s,abs(v(j)),beta(j),real(z(j)),imag(z(j))));
  end    
end
disp(' ')
if imag(c) < 0
  s = '-';
else
  s = '+';
end
disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))
disp(sprintf('\n  Conformal modulus = %.8g',imag(rect(2))/rect(1)/2));
%disp(sprintf('\n  Rectangle corners:'))
%R = [rect(1);real(rect(2));imag(rect(2));real(rect(3));imag(rect(3));rect(4)];
%disp(sprintf('    %.4f,  %.4f + %.4fi,  %.4f + %.4fi,  %.4f',R))

