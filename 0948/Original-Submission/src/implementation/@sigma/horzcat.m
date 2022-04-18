function c = horzcat(varargin)

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
    if row~=row2
        error('Cannot process object matrices of different rows!')
    end
    col = col + col2;
end
c(1:row, 1:col) = sigma(siz);

%% Assign values
colEnd = 0;
for i=1:n
    colStart = 1 + colEnd;
    colEnd = colEnd + size(varargin{i}, 2);
    if isa(varargin{i}, 'sigma')
        c(:, colStart:colEnd) = varargin{i};
    end
end
end
