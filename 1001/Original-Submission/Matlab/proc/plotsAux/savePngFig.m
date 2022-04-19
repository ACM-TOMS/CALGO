%% savePngFig
% Saves the figure as PNG, EPS or FIG.
%
%% Syntax
%
%   savePngFig(figNo,iOut,seti)
%
%% Description
% |savePngFig(figNo,iOut,seti)| saves the figure with number |figNo| at
% outer iteration |iOut| as PNG, EPS or FIG.
%
%% Example
%
%   figNo = 1;
%   iOut = 0;
%
%   seti.dirname = 'saveExample';
%   seti.fileSuffix = '_saveExample';
%
%   seti.savepng  = 1;
%   seti.saveepsc = 0;
%   seti.savefig  = 0;
%
%   figure(1);
%   x = -5:0.1:5;
%   plot(x,sin(x));
%
%   % Make dir saveExample (see above) if it does not exist
%   if exist(seti.dirname,'dir') ~= 7
%     mkdir(seti.dirname);
%   end
%
%   savePngFig(figNo,iOut,seti);
%
%% Input Arguments
%
% * figNo : figure number of the defined figure
% * iOut  : number of outer iteration (is used in filename)
% * seti  : structural array
%
% * seti.dirname      : dirname to save files, see <dirMake.html>.
% * seti.fileSuffix   : suffix in filename, see <setInput.html>.
%
% The following fields are also described in <setGeomSim.html> 
% in the section "Subfunction: setFigureSettings":
%
% * seti.savepng        :   (default: 1) 0 or 1, save figures as *.png.
% * seti.saveepsc       :   (default: 0) 0 or 1, save figures as colored *.eps.
% * seti.savefig        :   (default: 0) 0 or 1, save figures as *.fig.
%
%% Output Arguments
%
% Saved figures.
%
%% See Also
%
% * <dirMake.html>
% * <setInput.html>
% * <setGeomSim.html>
% * <savePngFigSimple.html>
%
%% Code

function savePngFig(figNo,iOut,seti)

filePng = sprintf('%s/fig_%02d_iOut_%02d%s.png',seti.dirname,figNo,iOut,seti.fileSuffix);
fileEpsc = sprintf('%s/fig_%02d_iOut_%02d%s.eps',seti.dirname,figNo,iOut,seti.fileSuffix);

warning('off','all');
if seti.savepng == 1
    figNostr = sprintf('-f%d',figNo);
    print(figNostr,'-dpng',filePng);
end

if seti.saveepsc == 1
    figNostr = sprintf('-f%d',figNo);
    print(figNostr,'-depsc',fileEpsc);
end

% Input filename with ending .png.
% Replace .png with .fig and store as .fig.

% e.g.:
% file = sprintf('%s/fig_9_iOut_%02d%s.png',seti.dirname,iOut,seti.fileSuffix);
% print('-f9','-dpng',file);
% file = sprintf('%s/fig_9_iOut_%02d%s.fig',seti.dirname,iOut,seti.fileSuffix);

if seti.savefig == 1
    fileFig = strrep(filePng,'.png','.fig');
    savefig(figNo,fileFig);
end
warning('on','all');

end
