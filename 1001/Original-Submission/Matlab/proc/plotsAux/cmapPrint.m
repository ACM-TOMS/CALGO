%% litman
% Define colormap for print.
%
%% Example
%
%   colormap(cmapPrint);
%
%% More About
%
% Default is the colormap litman. It can be changed to cmapPrint in the code.
%
%% See Also
%
% * <litman.html>
%
%% Code

function cmap = cmapPrint()

cmap = flipud(gray); % inverse of colormap 'gray'

end
