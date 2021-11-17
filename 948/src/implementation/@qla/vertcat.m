function c = vertcat(varargin)

n = length(varargin);
for i=1:n
    if isa(varargin{i}, 'qla')
        siz = varargin{i}(1).size;
        lenType = length(varargin{i}(1).type);
        break;
    end
end

[row, col] = size(varargin{1});
for i=2:n
    [row2, col2] = size(varargin{i});
    if col~=col2
        error('Cannot process object matrices of different columns!')
    end
    row = row + row2;
end
tmp = qla(siz); tmp = setType(tmp, zeros(1,lenType));
c(1:row, 1:col) = tmp;


rowEnd = 0;
for i=1:n
    rowStart = 1 + rowEnd;
    rowEnd = rowEnd + size(varargin{i}, 1);
    if isa(varargin{i}, 'qla')
        c(rowStart:rowEnd, :) = varargin{i};
    end
end
end