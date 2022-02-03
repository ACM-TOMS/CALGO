function A = subsasgn(A, S, B)

n = S.subs{:};
lenA = length(A);
lenB = length(B);
lenn = length(n);

if ~isa(A, 'sigma')
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
    if isa(B, 'sigma')
        siz = B.size;
        if lenB>1 && lenn~=lenB
            error('Size mismatched!');
        end
        if nMax>lenA
            start = min(lenA+1, nMax);
            A(start:nMax) = sigma(siz);
        end
        A(n) = B;
    else
        if nMax>lenA
            start = min(lenA+1, nMax);
            A(start:nMax) = sigma(A(1).size);
        end
    end
end
end