%% gridUp
% Up-scaling of a grid.
%
%% Syntax
%
%   A = gridUp(B,nROI,nInv)
%
%% Description
% |A = gridUp(B,nROI,nInv)| scales up the grid B of size nInv x nInv to A
% of size nROI x nROI. Upscaling for 3D-arrays is supported too.
%
%% Examples
%
% *2D: from 5 x 5 matrix to 8 x 8 matrix*
%
%   nInv = 5;
%
%   % B: test matrix (2D)
%   for i = 1:nInv
%       for j = 1:nInv
%           B(i,j) = i*10+j; % 2D
%       end
%   end
%
%   nROI = 8;
%
%   A = gridUp(B,nROI,nInv); % grid up-scaling: from nInv x nInv to nROI x nROI
%
% _Result:_
% From B up scaled to A.
% 
% B =
% 
%     11    12    13    14    15
%     21    22    23    24    25
%     31    32    33    34    35
%     41    42    43    44    45
%     51    52    53    54    55
% 
% A =
% 
%     11    12    12    13    14    14    15    15
%     21    22    22    23    24    24    25    25
%     21    22    22    23    24    24    25    25
%     31    32    32    33    34    34    35    35
%     41    42    42    43    44    44    45    45
%     41    42    42    43    44    44    45    45
%     51    52    52    53    54    54    55    55
%     51    52    52    53    54    54    55    55
%
%% Input Arguments
%
% * nInv    : length of matrix B
% * B in 2D : matrix of of size nInv x nInv
% * B in 3D : array of of size nInv x nInv x nInv
% * nROI    : length of resulting matrix A
%             (nROI x nROI in 2D, nROI x nROI x nROI in 3D)
%
%% Output Arguments
%
% * A in 2D: up scaled matrix size nROI x nROI
% * A in 3D: up scaled array of size nROI x nROI x nROI
%
%% More About
%
% The "inverse" of up-scaling is down-scaling function <gridDown.html>.
%
%% See Also
%
% * <gridDown.html>
%
%% Code

function A = gridUp(B,nROI,nInv)
% input B, output A
% grid up-scaling from nInv to nROI

if nROI < nInv
    A = B;
    disp('gridUp is not used because nROI must be equal or greater than nInv.')
else
    % for i = 1:nROI
    % 	c = ceil(i*nInv/nROI);
    % 	fprintf('i = %g | ceil = %g \n',i,c);
    % end

    dim = length(size(B));
    i = 1:nROI;
    c = ceil(i*nInv/nROI);

    %for j = 1:nROI
    %	A(1,j) = B(1,ceil(j*nInv/nROI));
    %end

    if dim == 2
        A = B(c,c);
    elseif dim == 3
        A = B(c,c,c);
    else
        error('Dimension has to be 2 or 3');
    end
    
end

end
