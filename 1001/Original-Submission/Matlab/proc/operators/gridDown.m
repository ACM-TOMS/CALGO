%% gridDown
% Down-scaling of a grid.
%
%% Syntax
%
%   B = gridDown(A,nROI,nInv)
%
%% Description
% |B = gridDown(A,nROI,nInv)| scales down the grid A of size nROI x nROI
% to B of size nInv x nInv. Downscaling for 3D-arrays is supported too.
%
%% Examples
%
% *2D: from 8 x 8 matrix to 5 x 5 matrix*
%
%   nROI = 8;
%
%   % A: test matrix (2D)
%   for i = 1:nROI
%       for j = 1:nROI
%           A(i,j) = i*10+j; % 2D
%       end
%   end
%
%   nInv = 5;
%
%   B = gridDown(A,nROI,nInv); % grid down-scaling: from nROI x nROI to nInv x nInv
%
% _Result:_
%
% From A down scaled to B.
% 
%   A =
% 
%     11    12    13    14    15    16    17    18
%     21    22    23    24    25    26    27    28
%     31    32    33    34    35    36    37    38
%     41    42    43    44    45    46    47    48
%     51    52    53    54    55    56    57    58
%     61    62    63    64    65    66    67    68
%     71    72    73    74    75    76    77    78
%     81    82    83    84    85    86    87    88
% 
% Downscaling with method 'imresize':
%
%   B =
% 
%     11    13    15    16    18
%     31    33    35    36    38
%     51    53    55    56    58
%     61    63    65    66    68
%     81    83    85    86    88
%
% Downscaling with method 'for':
%
%   B =
% 
%    11.0000   12.5000   14.0000   15.5000   17.5000
%    26.0000   27.5000   29.0000   30.5000   32.5000
%    41.0000   42.5000   44.0000   45.5000   47.5000
%    56.0000   57.5000   59.0000   60.5000   62.5000
%    76.0000   77.5000   79.0000   80.5000   82.5000
%
%% Input Arguments
%
% * nROI    : length of matrix A
% * A in 2D : matrix of size nROI x nROI
% * A in 3D : array of size nROI x nROI x nROI
% * nInv    : length of resulting matrix B
%             (nInv x nInv in 2D, nInv x nInv x nInv in 3D)
%
%% Output Arguments
%
% * B in 2D : down scaled matrix of of size nInv x nInv.
% * B in 3D : down scaled array of of size nInv x nInv x nInv.
%
%% More About
%
% The function |gridDown| uses the MATLAB function |imresize| with the
% interpolation method |nearest| in 2D case. The function |imresize| is
% applied on slices of the grid in case of 3D, see [1, Sec. 5.1].
%
% We have implemented two methods to scale down the grid, which have a 
% similar result, but differ strongly in computational time. We recommend 
% to use |method = 'imresize'| instead of |method = 'for'|.
%
% The "inverse" of down-scaling is up-scaling function <gridUp.html>.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <gridUp.html>
%
%% Code

function B = gridDown(A,nROI,nInv)

method = 'imresize'; % 'imresize': matlab routine; alternative: 'for': a slow for-loop

dim = length(size(A));

% test different methods in 3D
if dim == 3
    %method = 'for';
    method = 'imresize';
    % 'imresize' is faster than 'for' method, e.g. 0.046 s instead of 0.311 s
end

switch method

    %%
    % *Code: method "imresize" (fast)*
    case 'imresize'

        %numrows = nInv;
        %numcols = nInv;
        %B = imresize(A, [numrows numcols],'nearest');
        if dim == 2
            B = imresize(A, [nInv nInv],'nearest');
        elseif dim == 3
            C = zeros(nInv,nInv,nROI);
            for ii = 1:nROI
                C(:,:,ii) = imresize(A(:,:,ii), [nInv nInv],'nearest');
            end
            B = zeros(nInv,nInv,nInv);
            for ii = 1:nInv
                B(ii,:,:) = imresize(squeeze(C(ii,:,:)), [nInv nInv],'nearest');
            end
        else
            error('Dimension has to be 2 or 3');
        end

    %%
    % *Code: method "for" (slow)*
    case 'for' % 2D and 3D

        % grid down

        % slow in case of large matrices...

        %-- I down: input A, output B
        %grid down-scaling from nROI to nInv

        if nROI < nInv
            B = A;
            disp('gridDown is not used because nROI must be equal or greater than nInv.')
        else
            
            dim = length(size(A));
            i = 1:nROI;
            c = ceil(i*nInv/nROI);

            % [u,ic,iu] = unique(c); % iu = c
            [~,ic,~] = unique(c); % u and iu are not used

            % ic: index of begin of entries with number... (in octave end of...)

            % numbers of entries...:
            %e2 = circshift(ic',1)';
            %e2(1) = 0;
            %n = ic-e2; % number of elements from u
            % e.g.: c = 1 2 2 3 3; u = 1 2 3; n = 1 2 2; (counting: 1x 1; 2x 2; 2x 3)

            iS = zeros(1,nInv);
            iE = zeros(1,nInv);
            for ii = 1:nInv
                %if ii == 1; istart = 1; else; istart = ic(ii-1)+1; end;
                %alternative: extend e: first entry 1: and then iv = e(k):e(k+1)
                istart = ic(ii);
                if ii == nInv
                    iend = nROI;
                else
                    iend = ic(ii+1)-1;
                end
                %iv = istart:iend
                iS(ii) = istart;
                iE(ii) = iend;
            end

            % this for-loop is slow...
            if dim == 2
                B = zeros(nInv,nInv);
                for l = 1:nInv
                    i = iS(l):iE(l);
                    for m = 1:nInv
                        j = iS(m):iE(m);
                        B(l,m) = mean(mean(A(i,j)));
                    end
                end
            elseif dim == 3
                B = zeros(nInv,nInv,nInv);
                for l = 1:nInv
                    for m = 1:nInv
                        for n = 1:nInv
                            i = iS(l):iE(l); 
                            j = iS(m):iE(m);
                            k = iS(n):iE(n);
                            B(l,m,n) = mean(mean(mean(A(i,j,k)))); %in 3D: mean(mean(mean(...))), in 2D: mean(mean(...)) would be ok
                        end
                    end
                end
            else
                error('Dimension has to be 2 or 3');
            end
        end % if
end % switch

end
