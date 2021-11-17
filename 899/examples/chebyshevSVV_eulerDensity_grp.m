
% chebyshevSVV_eulerDensity_grp

       clear, home, close all

       alpha = 0.99;                                % KT map parameter
       % Note: depends on where the MPT is installed
       fN = fullfile('C:\docs\Matlab','suites','pp','PdeExamples','chebyshevEuler128.txt');   
       fa = dlmread(fN,'\n');  fa = fa(:)';                    
       N = length(fa);
       
        x = -cos((0:N-1)*pi/(N-1));       
       xa = asin(-alpha*cos(pi*(0:N-1)/(N-1)))/asin(alpha);   
       
       ak = chebyshevCoefficients(fa); 
       
       fN = fullfile('C:\docs\Matlab','suites','pp','PdeExamples','chebyshevEuler128exact.txt');
       fExact = dlmread(fN,'\n');  fExact = fExact(:)';
       
       S = [-1 -0.4456 -0.0222 0.3461 0.6606 1];       % mapped grid edge locations - discontinuities in f and f'
       S = sin( S*asin(alpha) )/alpha;                 % CGL grid edge locations to be passed to grp.m
       LK = [5 6 6 4 4];                               % Lamda in each smooth subinterval
       MK = [4 12 4 4 4];                              % m in each smooth subinterval

      [up,gh] = grp(S,LK,MK,ak,x,1);
      
      plot(xa,fa,'b',xa,up,'r')
      xlabel 'x', ylabel '\rho(x,t=0.4)'
      figure
      semilogy(xa,abs(up-fExact),'r')
      xlabel 'x', ylabel '|grp - exact|'
      
