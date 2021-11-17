% EXAMPLEFUNCTIONSETUP2D  Executes when a function is selected from the Examples2d menu of the MPT GUI.
%
% Functions called:
%  1) exampleFunctions2d.m
% Called by:
%  1) mpt.m
% Last modified: October 17, 2007

function exampleFunctionSetup2d(h, ha, fCh)         
 
     ha.isPdeSolution = false;
     ha.is2d = true;
     
     N = str2double(get(ha.N,'String'));
     M = str2double(get(ha.M,'String'));
   
    i=0:N;                                           
    polyCh = get(ha.chebyButton,'Value');
    if polyCh == 1                                     %  Chebyshev
       x = -cos(i*pi/N);
       xp = linspace(-1,1,M);
    else                                               %  Fourier
        x = -1 + 2*(0:N-1)/N;
       xp = -1 + 2*(0:M-1)/M;  
    end
     
       [X,Y] = meshgrid(x,x);
       [XP,YP] = meshgrid(xp,xp); 
   
    [f,fp] = exampleFunctions2d(X,Y,XP,YP,fCh);
    
       X = []; Y = []; XP = []; YP = [];
       
       switch polyCh
         case 1
           [uc,ak] = chebyshevInterpolation2d(f,xp,xp);
         otherwise
            [uc,ak] = fourierInterpolation2d(f,M,M);
       end
    
    if strcmp(get(ha.countouMenuItem, 'Checked'), 'on')     % axes1 - upper plot, axes2 - lower
      contour(ha.axes1,xp,xp,uc);
      rotate3d off
    else
      surf(ha.axes1,xp,xp,uc);
      rotate3d on
    end
  

    ha.xa = xp; 
    ha.f = f;
    ha.fExact = fp;
    ha.fa = uc;
    ha.ak = ak;
    ha.x = x;
    
    
    if strcmp(get(ha.countouMenuItem, 'Checked'), 'on')
      contour(ha.axes2,xp,xp,abs(uc-fp));
    else
      surf(ha.axes2,xp,xp,abs(uc-fp));
    end
    
    hlink = linkprop([ha.axes1 ha.axes2],'View');
    ha.hlink = hlink;
    guidata(h,ha);
