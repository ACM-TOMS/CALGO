function hvt = findHVT(A, r)

n_blocks = length(r) - 1;

n = size(A,1);
hvt = zeros(1,n);

for i = 1 : n_blocks
  % Get a diagonal sub-block
	rowcol = r(i):r(i+1)-1;
	B = A(rowcol, rowcol);
  
    if ~issparse(B)
        hvt(rowcol) = lapjv(double(int32(-B)));
    else
        hvt(rowcol) = lapjv(-B);
    end
end
end



