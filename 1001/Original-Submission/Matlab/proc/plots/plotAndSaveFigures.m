%% plotAndSaveFigures
% Plot and save figures of reconstruction.
%
% Reserved numbers: 11-20. Used numbers 11-16, see <start.html>.
%
%% Syntax
%
%   plotAndSaveFigures(seti,qROIexact,qROIcomp,iOut,out)
%
%% Description
% |plotAndSaveFigures(seti,qROIexact,qROIcomp,iOut,out)| plots and saves
% figures of outer iteration |iOut| in dependence of |out|.
%
% In case of 3D sectional planes through the contrast are plotted to. The
% sectional plane is prepared in the subfunction |plot3Dsec|.
%
%% Input Arguments
%
% * seti      : structural array
% * qROIexact : predefined (i.e. exact) contrast
% * qROIcomp  : reconstructed (i.e. computed) contrast
% * iOut      : number of outer iteration
% * out       : Output depth: generate no plots (0),
%               generate plots (1), generate plots and save them (2).
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
% The most important fields for this function in |seti| are:
%
% * seti.dis : relative discrepancies
% * seti.err : relative errors
% * seti.MTv : result of Tikhonov funtional
% * seti.M1v : result of first part of Tikhonov functional
% * seti.M2v : result of second part of Tikhonov functional
%
% Further details of this fields are described in <recon.html>.
%
%% Output Arguments
%
% Plots, see <start.html>.
%
% Note that figure 11 is plotted in <subplots.html>, but saved in this
% file.
%
%
%% See Also
%
% * <start.html>
% * <savePngFig.html>
% * <setGeomSim.html>
% * <recon.html>
% * <subplots.html>
% * <plot2DstylePublish.html>
% * <plot3DstylePublish.html>
% * <litman.html>
%
%% Code
%
function plotAndSaveFigures(seti,qROIexact,qROIcomp,iOut,out)
% qROIexact is seti.qROIexact.
% Note that this function will be used with source problem too
% (source problem is not available in public version).

%%
% *iOut > iOutIni: figure 11 (save subplots)*

if out == 2 && iOut > seti.iOutIni % iOut = iOutIni shows only the true contrast
    
    % figure 11: subplots (save)
    savePngFig(11,iOut,seti); % figure 11: subplots (was already plotted but not saved)
end

%%
% *iOut > iOutIni: figure 12-16*

if out >= 1 && iOut > seti.iOutIni

%%
% *figure 12: discrepancy and error (as in subplot)*

    figure(12);
    set(gcf,'Visible',seti.plotVisible);
    hplot = plot((0:iOut), [seti.disIni seti.dis(1:iOut)], 'b-', (0:iOut), [seti.errIni seti.err(1:iOut)], 'r:');
    if seti.plotPublish
        set(gca,'FontSize',seti.pubFontSize)
        lw = 2;
        set(hplot(1),'LineWidth',lw);
        set(hplot(2),'LineWidth',lw);
    else
        legend('discrepancy','error');
        title('Discrepancy and error');
    end
    axis square;
    if out == 2
        savePngFig(12,iOut,seti);
    end

%%
% *figure 13: Tikhonov functional (as in subplots)*

    figure(13);
    set(gcf,'Visible',seti.plotVisible);
    x = 1:iOut;
    if length(x) == 1
        [ax,h1,h2] = plotyy(x',seti.M1v(x)',x',seti.M2v(x));
    else
        [ax,h1,h2] = plotyy(x',seti.M1v(x)',[x',x'],[seti.M2v(x)',seti.MTv(x)']);
    end
    if seti.plotPublish
        set(gca,'FontSize',seti.pubFontSize)
        lw = 2;
        set(hplot(1),'LineWidth',lw);
        set(hplot(2),'LineWidth',lw);
    else
        if length(x) == 1
            legend([h1;h2],'M1 (left axis)','M2 (right axis)');
        else
            legend([h1;h2],'M1 (left ax.)','M2 (right ax.)','MT = M1 + M2 (right ax.)');
        end
        title('Min. Tikhonov funct. MT');
        xlabel('outer iterations iOut')
    end
    axis(ax,'square');
    if out == 2
        savePngFig(13,iOut,seti);
    end

%%
% *figure 14: reconstructed contrast (real part)*

    figure(14);
    set(gcf,'Visible',seti.plotVisible);
    if seti.dim == 2
        %imagesc(seti.G(real(qROIcomp))); colormap(cmapPrint); colorbar; % alternative colormap
        imagesc(seti.G(real(qROIcomp))); colormap(litman); colorbar;
        axis xy;
        if seti.usecbarlim
            caxis(seti.cbarlim);
        end
    elseif seti.dim == 3
        % isosurface(seti.G(real(qROIPrev))); colormap(litman); caxis([-1,3]); colorbar;
        % title('real part of reconstructed contrast (isosurface plot)');
        contourPlotROI(qROIcomp, seti, 'real');
    end
    if seti.plotPublish
        if seti.dim == 2
            plot2DstylePublish;
        elseif seti.dim == 3
            plot3DstylePublish;
        end
    else
        title('reconstructed contrast (real part)');
    end
    axis square;
    if out == 2
        savePngFig(14,iOut,seti);
    end
    close(14);
    % It is important to close this figure 
    % because otherwise the projection in 3D is a problem.

%%
% *figure 15: reconstructed contrast (imag part)*

    figure(15);
    set(gcf,'Visible',seti.plotVisible);
    if seti.dim == 2
        imagesc(seti.G(imag(qROIcomp))); colormap(litman); colorbar;
        axis xy;
    elseif seti.dim == 3
        contourPlotROI(qROIcomp, seti, 'imag');
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
        title('reconstructed contrast (imag part)');
    end
    axis square;
    if out == 2
        savePngFig(15,iOut,seti);
    end
    close(15);
    % It is important to close this figure 
    % because otherwise the projection in 3D is a problem.

%%
% *figure 16: difference of reconstructed and true contrast (abs)*

if 0
    figure(16);
    set(gcf,'Visible',seti.plotVisible);
    img = seti.G(abs(qROIprev-qROIexact));
    if seti.dim == 3
        % Display only a slice of the data:
        N2 = ceil((seti.nROI-1)/2);
        img = img(:,:,N2);
    end
    imagesc(img); colormap(litman); colorbar;
    axis xy;
    if seti.usecbarlim
        caxis(seti.cbarlim);
    end
    if seti.plotPublish
        plot2DstylePublish;
    else
        title('difference of reconstructed and true contrast (abs)');
    end
    axis square;
    if out == 2
        savePngFig(16,iOut,seti);
    end
end

end

%%
% *figures 21-30 reserved for 3D reconstruction*
%
% * figure 21-26
% * 21-23: predefined contrast (real part) (3D sections)
% * 24-26: reconstructed contrast (real part) (3D sections)
%
if out >= 1 && seti.dim == 3 && ...
   (strcmp('corner3D',seti.contrast) || strcmp('twoTripods3D',seti.contrast) ||...
    strcmp('twoTripodsRealImag3D',seti.contrast) || strcmp('cubeLike3D',seti.contrast) ||...
    strcmp('cross3D',seti.contrast))
    switch seti.contrast
        case {'corner3D', 'twoTripods3D', 'twoTripodsRealImag3D', 'cubeLike3D'}
            seca = -(1-3.5/8);
            %secn = 2; % most interesting
        case 'cross3D'
            seca = 0;
            %secn = 2; % 1, 2 or 3... in coordinate system
    end
    if iOut == 0
        plot3Dsec(qROIexact,'predefined',seca,1,seti,iOut,out,21);
        plot3Dsec(qROIexact,'predefined',seca,2,seti,iOut,out,22);
        plot3Dsec(qROIexact,'predefined',seca,3,seti,iOut,out,23);
    else
        plot3Dsec(qROIcomp,'reconstructed',seca,1,seti,iOut,out,24);
        plot3Dsec(qROIcomp,'reconstructed',seca,2,seti,iOut,out,25);
        plot3Dsec(qROIcomp,'reconstructed',seca,3,seti,iOut,out,26);
    end
end

end

%%
% *Code: subfunction: plot3Dsec*
%
% plot sectional plane in case of 3D

function plot3Dsec(qROI,contrastKind,seca,secn,seti,iOut,out,figno)
% qROI: seti.qROIexact (exact/predefined) OR qROIprev (reconstructed)
% contrastKind: in title of plot, i.e. predefined or reconstructed
R = seti.rCD/2;
ind = floor((seca*R+R)/(2*R)*seti.nROI);
A = seti.G(qROI);

switch secn
    case 1
        B = A(ind,:,:);
    case 2
        B = A(:,ind,:);
    case 3
        B = A(:,:,ind);
    otherwise
        error('n must be 1, 2 or 3.');
end
B = squeeze(real(B));
B = transpose(B);

figure(figno);
set(gcf,'Visible',seti.plotVisible);
imagesc(B); colormap(litman); colorbar;
axis xy;
if seti.usecbarlim
    caxis(seti.cbarlim);
end
if seti.plotPublish
    plot2DstylePublish;
else
    titles = sprintf('%s contrast (real part), X%d',contrastKind,secn);
    title(titles);
end
axis square;
if out == 2
    savePngFig(figno,iOut,seti);
end
end
