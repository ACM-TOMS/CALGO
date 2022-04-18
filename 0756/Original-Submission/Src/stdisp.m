function stdisp(w,beta,z,c)
%STDISP Display results of Schwarz-Christoffel strip parameter problem.
%       STDISP(W,BETA,Z,C) displays the results of STPARAM in a pleasant
%       way.  
%
%	See also STPARAM, STPLOT.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

disp(' ')
disp('          w                beta                  z           ')
disp(' ------------------------------------------------------------')
u = real(w);
v = imag(w);
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  if ~imag(z(j))
    disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
      u(j),s,abs(v(j)),beta(j),z(j)));
  else
    disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e + i',...
      u(j),s,abs(v(j)),beta(j),z(j)));
  end    
end
disp(' ')
if imag(c) < 0
  s = '-';
else
  s = '+';
end
disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))

