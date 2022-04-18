function A = subsasgn(A, S, B)

dim = length(S.subs);
global DAEZERO;

if dim==1
    n = S.subs{:};
    lenA = length(A);
    lenB = length(B);
    lenn = length(n);
    
    if ~isa(A, 'qla')
        clear A;
    end
    
    nMax = max(n);
    
    if lenB==0
        if nMax>lenA, error('Subscript exceeds bound!'); end
        temp = A;
        clear A;
        diffset = setdiff(1:lenA, n);
        A(1:(lenA-lenn)) = temp(diffset);
    else
        if isa(B, 'qla')
        siz = B.size;
        if lenB>1 && lenn~=lenB
            error('Size mismatched!');
        end
        if nMax>lenA
            start = min(lenA+1, nMax);
            A(start:nMax) = qla(siz);
        end
        A(n) = B;
        
        else
            A(n) = DAEZERO;
        end
    end
    
elseif dim==2
    rows = S.subs{1};
    cols = S.subs{2};
    numrow = length(rows);
    numcol = length(cols);
    siz = B(1).size;
    
    if ~isa(A, 'qla')
        clear A;
        A(1:rows(end),1:cols(end)) = qla(siz);
    else
        [m, n] = size(A);
        A(min(m, rows(1)):rows(end), :) = qla(siz);
        A(:, min(n, cols(1)):cols(end)) = qla(siz);
    end
    
    if numel(B)>1
        [m, n] = size(B);
        if numrow~=m || numcol~=n
            error('Size mismatched!');
        end
        for i=1:numrow
            for j=1:numcol
                A(rows(i),cols(j)) = B(i,j);
            end
        end
    else
        for i=1:numrow
            for j=1:numcol
                A(rows(i),cols(j)) = B;
            end
        end
    end
else
    error('Dimension not supported!');
end
end