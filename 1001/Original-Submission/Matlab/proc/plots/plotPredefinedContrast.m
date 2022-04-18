%% plotPredefinedContrast
% Plot and save predefined contrast (real and imaginary part) in figure 2
% and 3.
%
%% Syntax
%
%   plotPredefinedContrast(seti,out)
%
%% Description
%
% |plotPredefinedContrast(seti,out)| is called in <setContrast.html> to plot
% and save (in case of |out = 2|) the real and imaginary part of the predefined contrast
% |seti.qROIexact| in figure 2 and 3.
%
%% Input Arguments
%
% * seti    :   structural array
% * out     :   output depth: no figure (0), plot figure (1), 
%               plot and save figure (2).
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
%% Output Arguments
%
% * figure 02: predefined contrast (real part)
% * figure 03: predefined contrast (imag part)
%
%% More About
%
% Note that |seti.G| is used to write the vector as a matrix in ROI.
%
%% See Also
% * <start.html>
% * <setContrast.html>
% * <contourPlotROI.html>
% * <setGeomSim.html>
% * <savePngFig.html>
% * <litman.html>
%
%% Code
%
function plotPredefinedContrast(seti,out)

%%
% *figure 2: predefined contrast (real part)*

figure(2);
set(gcf,'Visible',seti.plotVisible);
if seti.dim == 2
    imagesc(seti.G(real(seti.qROIexact))); colormap(litman); colorbar;
    axis xy;
elseif seti.dim == 3
    contourPlotROI(seti.qROIexact, seti, 'real');
end
if seti.usecbarlim
    caxis(seti.cbarlim);
end
if seti.plotPublish
    if seti.dim == 2
        plot2DstylePublish;
    elseif seti.dim == 3
        plot3DstylePublish;
    end
else
    title('predefined contrast (real part)');
end
axis square;
if out == 2
    savePngFig(2,0,seti); % set iOut = 0
end

%%
% *figure 3: predefined contrast (imag part)*

figure(3);
set(gcf,'Visible',seti.plotVisible);
if seti.dim == 2
    imagesc(seti.G(imag(seti.qROIexact))); colormap(litman); colorbar;
    axis xy;
elseif seti.dim == 3
    contourPlotROI(seti.qROIexact, seti, 'imag');
end
if seti.usecbarlim
    caxis(seti.cbarlim);
end
if seti.plotPublish
    if seti.dim == 2
        plot2DstylePublish;
    elseif seti.dim == 3
        plot3DstylePublish;
    end
else
    title('predefined contrast (imag part)');
end
axis square;
if out == 2
    savePngFig(3,0,seti); % set iOut = 0
end

end
