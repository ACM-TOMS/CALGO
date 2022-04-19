%% extendROIto CD
% Extends ROI to a CD.
%
%% Syntax
%
%   y = extendROItoCD(x, ROImask)
%
%% Description
% |y = extendROItoCD(x, ROImask)| extends the vector (or matrix) |x| in 
% the region of interest (ROI) to the computational domain (CD) by zeros, 
% where |ROImask| is logical to describe ROI inside CD.
%
%% Examples
%
% *Example 1: x is a vector*
%
%   ROImask = [0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0];
%   x = [1 2 3 4];
%   extendROItoCD(x,ROImask)                       
% 
% _Result:_
%
%   ans =
% 
%      0     0     0     0
%      0     1     3     0
%      0     2     4     0
%      0     0     0     0
%
% Note that this is another result than in Example 2.
%
% *Example 2: x is a matrix*
%
%   ROImask = [0 0 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0];
%   x = [1 2; 3 4];
%   extendROItoCD(x,ROImask)
% 
% _Result:_
%
%   ans =
% 
%      0     0     0     0
%      0     1     2     0
%      0     3     4     0
%      0     0     0     0
%
%% Input Arguments
% 
% * x   :   a vector or a matrix, in this framework
%           x is the contrast (e.g. seti.qROIexact) on ROI 
%           written as vector of size nROI^seti.dim x 1.
% * ROImask :   mask (logical matrix of size seti.nCD x seti.nCD)
%               to describe ROI inside CD by 1-entries; 
%               in this framework we define seti.ROImask, see <setGrid.html>
%
%% Output Arguments
%
% * y : ROI extended to CD by zeros.
%
%% See Also
%
% * <restrictCDtoROI.html>
%
%% Code
function y = extendROItoCD(x, ROImask)
y = zeros(size(ROImask)); % (size of CD is computed from ROImask)
y(ROImask~=0) = x;
end
