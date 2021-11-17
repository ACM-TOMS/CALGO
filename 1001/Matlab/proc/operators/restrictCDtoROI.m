%% restrictCDtoROI
% Restricts the CD to ROI.
%
%% Syntax
%
%   y = restrictCDtoROI(x,ROImask)
%
%% Description
% |y = restrictCDtoROI(x, ROImask)| restricts the vector (or matrix) |x| in 
% the computational domain (CD) to the region of interest (ROI), 
% where |ROImask| is logical to describe ROI inside CD.
%
%% Examples
%
% *Example 1: x is a vector*
%
%   x = [0 0 0 0 0 1 2 0 0 3 4 0 0 0 0 0];
%   ROImask = logical([0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]);
%   y = restrictCDtoROI(x,ROImask);
%
% _Result:_
%
%   y =
%
%        1     2     3     4
%
% *Example 2: x is a matrix*
%
%   x = [0 0 0 0; 0 1 2 0; 0 3 4 0; 0 0 0 0];
%   ROImask = logical([0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0]);
%   y = restrictCDtoROI(x,ROImask);
%
% _Result:_
%
%   y =
%
%     1
%     3
%     2
%     4
%
%
%% Input Arguments
% 
% * x       :   a vector or a matrix; in this framework 
%               x is the contrast in CD 
%               written as vector of size nCD^seti.dim x 1.
% * ROImask :   mask (logical matrix of size seti.nCD x seti.nCD) 
%               to describe ROI inside CD by 1-entries; 
%               in this framework we define seti.ROImask, see <setGrid.html>
%
%% Output Arguments
%
% * y : CD restricted to ROI 
%       (The resut is always a vector of size nROI^seti.dim x 1.)
%
%% See Also
%
% * <extendROItoCD.html>
%
%% Code
function y = restrictCDtoROI(x, ROImask)
y = x(ROImask(:));
end
