% SETUPPDEEXAMPLE  Executes when a selcted is made from the PDE Examples menu of the MPT GUI.
%
% Called by:
%  1) mpt.m
%
% Last modified: July 9, 2008
 
function setUpPdeExample(h,ha,gCh,alpha,fileName,fileNameExact)   

       ha.isPdeSolution = true;
       ha.alpha = alpha;
        
       fN = fileName;
       ha.fa = dlmread(fN,'\n');  ha.fa = ha.fa(:)';                    % ensure it is a row vector
       N = length(ha.fa);

      set(ha.N,'String',num2str(N));
      set(ha.M,'String',num2str(N));
      
      ha.f = ha.fa;
   
      polyCh = get(ha.chebyButton,'Value'); 
      switch polyCh
            case 1
                ha.x = -cos((0:N-1)*pi/(N-1));       % CGL grid
                ha.xa = asin(-alpha*cos(pi*(0:N-1)/(N-1)))/asin(alpha);   % CGL grid or a mapping of the CGL grid
                ha.ak = chebyshevCoefficients(ha.f); 
            otherwise
                ha.x = -1 + 2*(0:N-1)/N;
                ha.xa = -1 + 2*(0:N-1)/N; 
                ha.ak = fft(ha.f);
        end
      
      fN = fileNameExact;
      ha.fExact = dlmread(fN,'\n');  ha.fExact = ha.fExact(:)';
   
      axes(ha.axes1)
      plot(ha.xa,ha.fExact,'r',ha.xa,ha.fa,'g')
      if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('exact','numerical'); end
      
      axes(ha.axes2)                                                 % plot error error in botton image on semilogy plot
      semilogy(ha.xa(2:end),abs(ha.fExact(2:end) - ha.fa(2:end)),'g')
      if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('error numerical'); end
      
      guidata(h,ha);

  
