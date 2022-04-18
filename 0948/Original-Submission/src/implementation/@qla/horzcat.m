function c = horzcat(varargin)

n = length(varargin);
%% Obtain size of qla
for i=1:n
    if isa(varargin{i}, 'qla')
        siz = varargin{i}(1).size;
        lenType = length(varargin{i}(1).type);
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
tmp = qla(siz); tmp = setType(tmp, zeros(1,lenType));
c(1:row, 1:col) = tmp;

%% Assign values
colEnd = 0;
for i=1:n
    colStart = 1 + colEnd;
    colEnd = colEnd + size(varargin{i}, 2);
    if isa(varargin{i}, 'qla')
        c(:, colStart:colEnd) = varargin{i};
    end
end
end
