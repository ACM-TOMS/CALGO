function hpdisp(w,beta,x,c)
%HPDISP Display results of Schwarz-Christoffel half-plane parameter problem.
%       HPDISP(W,BETA,X,C) displays the results of HPPARAM in a pleasant
%       way.  
%
%	See also HPPARAM, HPPLOT.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

if length(x) < length(w)
  x = [x(:);Inf];
end
disp(' ')
disp('          w                beta               x          ')
disp(' --------------------------------------------------------')
u = real(w);
v = imag(w);
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
      u(j),s,abs(v(j)),beta(j),x(j)));
end
disp(' ')
if imag(c) < 0
  s = '-';
else
  s = '+';
end
disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))

