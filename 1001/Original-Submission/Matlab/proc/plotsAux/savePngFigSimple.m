%% savePngFigSimple
% Saves a figure as PNG.
%
%% Syntax
%
%   savePngFig(figNo,dirname,filename,number)
%
%% Description
% |savePngFigSimple(figNo,dirname,filename,number)|
% saves the figure with number |figNo| in the directory |dirname| as PNG
% with filename composited of |number| and |filename|.
%
%% Example
%
%   figNo = 1;
%   dirname = 'saveExample';
%   filename = 'saveExample';
%   number = 0;
%
%   figure(1);
%   x = -5:0.1:5;
%   plot(x,sin(x));
%
%   % Make dir saveExample (see above) if it does not exist
%   if exist(dirname,'dir') ~= 7
%     mkdir(dirname);
%   end
%
%   savePngFigSimple(figNo,dirname,filename,number);
%
% _Result:_
%
% Creates a directory saveExample and saves the figure 1 in the file 
% |00000_saveExample.png|.
%
%% Input Arguments
%
% * figNo    : number of figure
% * dirname  : name of directory to save the file
% * filename : filename consits of this
% * number   : is a part of the resulting filename
%
%% Output Arguments
%
% Saved figures.
%
%% See Also
%
% <savePngFig.html>
%
%% Code
function savePngFigSimple(figNo,dirname,filename,number)

%filePng = sprintf('figiPda/%s/%s_%02d_%s.png',dirname,dirname,number,filename);
filePng = sprintf('%s/%05d_%s.png',dirname,number,filename);

figNostr = sprintf('-f%d',figNo);
print(figNostr,'-dpng',filePng);

end
