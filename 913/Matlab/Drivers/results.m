%
% Show output for IDRS testproblems, 
% Called by idrs_ex*.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

c = ['b','g','r','c','m','k','y'];
colour = colour+1;
if exist('flag')
   if ( flag == 1 )
      disp('Flag = 1, maximum number of iterations reached!');
   elseif ( flag == 2 )
      disp('Flag = 2, true relative residual above tolerance!');
   elseif ( flag == 3 )
      disp('Flag = 3, Break down!');
   end 
end
if exist('relres')
   disp(['|b - Ax|/|b| = ', num2str(relres)]);
end
if exist('iter')
   disp(['Iterations: ',num2str(iter)]);
end
if exist('resvec')
   xas = [0:1:length(resvec)-1];
   yval = log10(resvec/resvec(1));
   plot(xas,yval,c(colour));
end
if exist('replacements')
   disp(['Residual replacements: ',num2str(replacements)]);
end
disp(' ');
