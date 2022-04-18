% EXAMPLEFUNCTIONSETUP1D  Executes when a function is selected from the Examples1d menu of the MPT GUI.
%
% Functions called:
%  1) exampleFunctions1d.m
% Called by:
%  1) mpt.m
% Last modified: October 17, 2007

function exampleFunctionSetup1d(h, ha, fCh)        


    ha.is2d = false;
    ha.isPdeSolution = false;
    rotate3d off

     N = str2double(get(ha.N,'String'));
     M = str2double(get(ha.M,'String'));


                                                %  Functions are known at N+1 CGL pts (Cheby approx) or uniform (Fourier approx)
    if get(ha.chebyButton,'Value') == 1         %  Chebyshev
        x = -cos((0:N)*pi/N);
       xp = linspace(-1,1,M);
    else                                        %  Fourier
       x = -1 + 2*(0:N-1)/N;
       xp = -1 + 2*(0:M-1)/M;
    end

    [f,fp,fh] = exampleFunctions1d(x,xp,fCh);

    if get(ha.chebyButton,'Value') == 1               %  Chebyshev
        [uc,ak] = chebyshevInterpolation(f,xp);
    else                                               %  Fourier
      [uc,ak] = fourierInterpolation(f,xp);
    end

    axes(ha.axes1)                                      % plot exact and approximation in upper image
    plot(xp,fp,'r',xp,uc,'g')
    if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('exact','numerical'); end

    ha.xa = xp;
    ha.f = f;
    ha.fExact = fp;
    ha.fa = uc;
    ha.x = x;
    ha.ak = ak;
    ha.exampleHandle = fh;
    guidata(h,ha);

    axes(ha.axes2)                                                 % plot interpolation error in botton image on semilogy plot

    semilogy(ha.xa(2:end),abs(ha.fExact(2:end) - ha.fa(2:end)),'g');
    if strcmp(get(ha.displayLegends, 'Checked'), 'on'), legend('error numerical'); end
