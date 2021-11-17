function c = vertcat(varargin)

n = length(varargin);
%% Obtain size of sigma
for i=1:n
    if isa(varargin{i}, 'sigma')
        siz = varargin{i}(1).size;
        break;
    end
end

%% Determine size of output
[row, col] = size(varargin{1});
for i=2:n
    [row2, col2] = size(varargin{i});
    if col~=col2
        error('Cannot process object matrices of different columns!')
    end
    row = row + row2;
end
c(1:row, 1:col) = sigma(siz);

%% Assign values
rowEnd = 0;
for i=1:n
    rowStart = 1 + rowEnd;
    rowEnd = rowEnd + size(varargin{i}, 1);
    if isa(varargin{i}, 'sigma')
        c(rowStart:rowEnd, :) = varargin{i};
    end
end
end