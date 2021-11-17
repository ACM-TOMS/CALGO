%% litman
% Define colormap litman.
%
%% Example
%
%   colormap(litman);
%
%% See Also
%
% * <cmapPrint.html>
%
%% Code

function cmap = litman()

N = 256;
y = (0:N-1)'/(N-1);

% Original
R = interp1([.00 .25 .44 .57 .76 1.0], [.00 .00 .73 1.0 1.0 1.0],y);
G = interp1([.00 .25 .44 .57 .76 1.0], [.00 .00 .00 .25 1.0 1.0],y);
B = interp1([.00 .25 .44 .57 .76 1.0], [.00 1.0 1.0 .85 .42 1.0],y);

% Smoothed
%R = interp1([.00 .25 .50 .75 1.0], [.00 .00 1.0 1.0 1.0],y);
%G = interp1([.00 .25 .50 .75 1.0], [.00 .00 .25 1.0 1.0],y);
%B = interp1([.00 .25 .50 .75 1.0], [.00 1.0 .75 .25 1.0],y);

cmap = [R G B];

end
