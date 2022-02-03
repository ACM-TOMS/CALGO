%	  * Schwarz-Christoffel Toolbox demonstrations *
%
%		1)  Tutorial
%		2)  Infinite vertices
%               3)  Elongated polygons
%		4)  Faber polynomials
%
%		0)  Quit
%
%    Warning: All current workspace variables will be lost.
echo off

while 1
  demos = str2mat('tutdemo','infdemo','elongdemo','faberdemo');
  clc
  help scdemo
  n = input('Select a demo number: ');
  if ((n <= 0) | (n > 4)) 
    break
  end
  eval(demos(n,:))
  clear
end


