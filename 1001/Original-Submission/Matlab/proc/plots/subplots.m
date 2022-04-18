%% subplots
% Plotting subplots in figure 11: 
% (1) true contrast, (2) reconstructed contrast, (3) discrepancy and error,
% (4) Tikhonov functional minimization.
%
%% Syntax
%
%   subplots(seti,qROIexact,qROIcomp,iOut)
%
%% Description
% |subplots(seti,qROIexact,qROIcomp,iOut)| plots the four subplots in
% figure 11 for outer iteration |iOut|:
%
% * Subplot 1: true contrast
% * Subplot 2: reconstructed contrast
% * Subplot 3: discrepancy and error
% * Subplot 4: Tikhonov functional minimization
%
% Note that the command to save figure 11 is in <plotAndSaveFigures.html>.
%
%% Input Arguments
%
% * seti      : structural array
% * qROIexact : predefined (i.e. exact) contrast
%               (usually qROIexact = seti.qROIexact, but in other cases 
%               (e.g. 'source') you can put in something else to plot.)
% * qROIcomp  : reconstructed (i.e. computed) contrast
% * iOut      : number of outer iteration
%
% * seti.recname : reconstruction name, e.g. 'contrast' or 'source'.
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
%% Output Arguments
%
% Output is figure 11 as described in description and in section 
% "Output: Figures" in <start.html>.
%
%% See Also
%
% * <plotAndSaveFigures.html>
% * <contourPlotROI.html>
% * <setGeomSim.html>
% * <litman.html>
%
%% Code
function subplots(seti,qROIexact,qROIcomp,iOut)

figureSupplots = figure(11);
set(gcf,'Visible',seti.plotVisible);
clf;
sprows = 2; % subplot rows
spcols = 2; % subplot cols

%%
% *Subplot 1: true contrast*

subplot(sprows,spcols,1,'Parent',figureSupplots);
gridax = seti.gridROI(2,1:seti.nROI);
switch seti.dim
    case 2
        imagesc(gridax,gridax,seti.G(real(qROIexact))); colormap(litman); colorbar;
        axis xy;
        if seti.usecbarlim
            caxis(seti.cbarlim);
        end
        % imagesc(gridax,gridax,seti.G(real(seti.wW(seti.qROIexact)))); colormap(litman); colorbar;
        % to control if matrix of wavelet-coeff. is sparse
    case 3
        contourPlotROI(qROIexact, seti, 'real');
end
axis square;
title(sprintf('True %s (real part) on ROI \n (nROI = %g) (rel. noise \\delta = %g)',seti.recname,seti.nROI,seti.delta));

if iOut == 0
    return % other subplots are empty...
end

%%
% *Subplot 2: contrast reconstruction*

subplot(sprows,spcols,2,'Parent',figureSupplots);

        if seti.dim == 2
            gridax = seti.gridROI(2,1:seti.nROI);
            qROIplot = qROIcomp;
            imagesc(gridax,gridax,seti.G(real(qROIplot))); colormap(litman); colorbar;
            axis xy;
            if seti.usecbarlim
                caxis(seti.cbarlim);
            end
        elseif seti.dim == 3
            % gridax = seti.gridROI(2,1:seti.nROI);
            % isosurface(gridax,gridax,gridax,seti.G(real(qROI)));
            contourPlotROI(qROIcomp, seti, 'real');
        else
            disp('seti.dim = 2 or 3...')
        end
        title(sprintf('Reconstructed %s \n (real part)',seti.recname));
        axis square;

%%
% *Subplot 3: discrepancy and error*

subplot(sprows,spcols,3,'Parent',figureSupplots);
        plot((0:iOut), [seti.disIni seti.dis(1:iOut)], 'b-', (0:iOut), [seti.errIni seti.err(1:iOut)], 'r:');
        legend('discrepancy','error');
        title('Discrepancy and error');
        xlabel('outer iterations iOut')
        axis square;

%%
% *Subplot 4: minimize functional*

subplot(sprows,spcols,4,'Parent',figureSupplots);
    x = 1:iOut;
    if length(x) == 1
        [ax,h1,h2] = plotyy(x',seti.M1v(x)',x',seti.M2v(x));
        legend([h1;h2],'M1 (left axis)','M2 (right axis)');
    else
        [ax,h1,h2] = plotyy(x',seti.M1v(x)',[x',x'],[seti.M2v(x)',seti.MTv(x)']);
        legend([h1;h2],'M1 (left ax.)','M2 (right ax.)','MT = M1 + M2 (right ax.)');
    end
    title('Min. Tikhonov funct. MT');
    axis(ax,'square');
    xlabel('outer iterations iOut')
 
    % plot((1:iOut), seti.MTv(1:iOut), (1:iOut), seti.M1v(1:iOut), (1:iOut), seti.M2v(1:iOut));
    % legend('MT = M1 + M2','M1','M2');
    
    %%
    % Error in legend in case of first outer step, i.e. iOut == 1:
    
    % Same with this example code:
    % x = 1;
    % A(x) = 1;
    % B(x) = 2;
    % C(x) = 3;
    % [ax,h1,h2] = plotyy(x',A(x)',[x',x'],[B(x)',C(x)']);
    % legend([h1;h2],'A','B','C');
    %
    % Warning: Ignoring extra legend entries. 
    % > In legend>set_children_and_strings (line 643)
    %   In legend>make_legend (line 328)
    %   In legend (line 257) 
    %
    % No problem if x is a vector: e.g. x = 1:2;
    %
    % Kind of solution:
    % Do not plot MT for first outer iteration...

%%

% warning off to avoid: "Warning: Error updating Light. Exceeded the maximum number (8) of light sources"
warning off;
drawnow;
warning on;

end
