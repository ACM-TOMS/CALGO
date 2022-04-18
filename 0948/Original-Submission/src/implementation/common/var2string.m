function s = var2string(varsym, A, breaknum, spaceNum)

if nargin<=3, spaceNum=0; maxLen = []; end

[m, n] = size(A);
i = 1;

if ~isa(varsym, 'cell')
    % a cell array of length n
    sv = strcat(varsym, num2str(A(i,1)), order2string(A(i,2)));    
    s  = sprintf('%s', sv);    
    for i=2:m
        % create variable name plus derivative
        sv = strcat(varsym, num2str(A(i,1)), order2string(A(i,2)));
        s  = sprintf('%s, %s', s, sv);
    end
    
else
    % a single string or char
    sv = strcat(varsym{A(i,1)}, order2string(A(i,2)));
    s  = sprintf('%s', sv);    
    for i=2:m
        sv = strcat(varsym{A(i,1)}, order2string(A(i,2)));
        s  = sprintf('%s, %s', s, sv);
    end    
end

if nargin>=3
    s = breakLine(s, breaknum, ' ' ,spaceNum);
    % delimiter could also be ','
end
end