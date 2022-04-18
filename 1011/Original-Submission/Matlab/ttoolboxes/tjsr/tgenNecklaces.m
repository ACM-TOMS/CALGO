function [ X ] = tgenNecklaces( n, k)
% [ X ] = tgenNecklaces( n, k )
% tgenNecklaces generates all short n-bead necklaces with k colors
%    i.e., all equivalence classes of binary n-tuples under 
%           1) cyclic permutation (rotation) and
%           2) powers.
%    Each equivalence class being represented by the lexicographically first tuple in the class.
%
% Input:
%   n     length of necklace
%   k     number of colours
%
% Output:
%   X     n x L vector, all short n-bead necklaces
%
% E.g.: tgenNecklaces( 3, 3 )

% Original version: October 2018
% Modified by: tommsch, February 2019
%              tommsch, March 2019           Behaviour change


% Generate all k-ary N-tuples
A = zeros( k^n, n );
for idx = 1:n;
    ww = k^(n-idx);
    for i = 0:(k^idx)-1;
        A(i*ww+1:(i+1)*ww, idx) = mod(i, k); end; end;

% Remove all non-necklaces
scaler = k.^(n-1:-1:0)';
ok = true( k^n, 1 );
shifter = zeros( n-1, n );
for i = 1:n-1;
    shifter(i,:) = [(i+1):n, 1:i]; end;
for scan = 1:k^n;
    if( ok(scan) );
        base = A(scan, :);
        ok(1+base(shifter)*scaler) = 0;
        ok(scan) = 1; end; end;
X = A(ok,:).' + 1;

idx = true( 1, size(X,2) );
for i = 1:size(X,2);
    if( ~isequal(reducelength(X(:,i)),X(:,i)) ); 
        idx(i) = false; end; 
end;
X = X(:,idx);

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   
