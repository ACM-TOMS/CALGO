function x = kronm(Q,x)
% Fast Kronecker matrix multiplication, for both full and sparse matrices 
% of any size. Never computes the actual Kronecker matrix and omits
% multiplication by identity matrices.
% y = kronm(Q,x) computes
%     y = (Q{k} kron ... Q{2} kron Q{1})*x
% If Q contains only two matrices and x is a vector, the code uses the
% identity
%     ( Q{2} kron Q{1} )*vec(X) = vec(Q{1}*X*Q{2}'),
% where vec(X)=x. If Q contains more than two matrices and/or if x has more
% than one column, the algorithm uses a generalized form of this identity.
% The idea of the algorithm is to see x as a multi-dimensional array and to
% apply the linear maps Q{i} separately for each dimension i. If Q contains
% just one matrix, the function returns the regular matrix product Q{1}*x.
%
% Inputs:
% Q:        1-by-k cell array containing k matrices of arbitrary size 
%           (can be sparse). Denote by R(i) the number of rows of Q{i}, and 
%           by C(i) the number of columns. Alternatively, Q{i} may also be
%           a scalar qi. This is interpreted as the qi-by-qi identity
%           matrix. Hand over identity matrices in this fashion for optimal
%           performance.
% x:        Matrix of size CC-by-m, where CC=C(1)*...*C(k).
%
% Output:   Matrix of size RR-by-m, where RR=R(1)*...*R(k).
%
%
% Example:
% R = [60, 30, 20];           % Number of rows for matrices Q{1},Q{2},Q{3}.
% C = [55, 25, 15];           % Number of columns of matrices Q{i}.
% m = 5;                      % Number of columns of x.
% Q = cell(1,length(R));      % Create cell with sparse random matrices
% for i=1:length(R)           % of density 0.05.
%     Q{i} = sprand(R(i),C(i),0.05);
% end
% x = rand(prod(C),m);        % Random matrix x with C(1)*C(2)*C(3) rows.
% y = kron(Q{3},kron(Q{2},Q{1}))*x;
%                             % Matlab's Kronecker multiplication...
% yy= kronm(Q,x);             % and kronm...
% norm(y-yy)                  % ... give the same result up to 
%                             % computational error.
%
%
% Version: 6-Oct-2015
% Author:  Matthias Kredler (Universidad Carlos III de Madrid)
%          mkredler@eco.uc3m.es
% Acknowledgement:
% This code follows the same idea as 'kronmult' by Paul G. Constantine & 
% David F. Gleich (Stanford, 2009). However, I avoid loops and allow for
% non-square inputs Q{i}. I have also included the special treatment for
% identity matrices. 

m = size(x,2);                      % Obtain number of columns in input.
k = length(Q);                      % Number of matrices in Q.
R = zeros(1,k);                     % Vector for number of rows of,
C = zeros(1,k);                     % Q-matrices and for number of columns.
comp = true(1,k);                   % Check if we have to multiply by Q{i}.
for i=1:k
    if isscalar(Q{i})               % If input Q{i} is a scalar, don't 
       comp(i) = false;             % have to multiply in this dimension.
       R(i) = Q{i};                 % Read  in number of rows and columns.
       C(i) = Q{i};
    else                            % Otherwise, read out size of the 
       [R(i),C(i)] = size(Q{i});    % matrix.
    end
end

xsiz = [C,m];                       % Will constantly change dimension of x. 
                                    % xsiz is the current size, when x is
                                    % reshaped to array of dim.
                                    % C(1),C(2),...,C(k),m.

if comp(1)                          % Start with first Kronecker product,                               
    x = Q{1}*reshape(x,[C(1),prod(xsiz)/C(1)]);  
                                    % leave out if Q{i} is identity.
    xsiz(1) = R(1);                 % Replace size of dimension 1.
    
           % disp(size(Q{1}))
end                                 % (Don't do this in loop below --> save
                                    % time on reshapes and permutes)
if k>1 && m==1                      % If Q has just one element, we're done.
    if comp(k)                      % If x was a column vector, do the last
        x = reshape(x,[prod(xsiz)/C(k),C(k)]) *Q{k}' ;
        xsiz(k) = R(k);             % Kronecker product by matrix
        
            %disp(size(Q{k}))
    end                             % post-multiplication to save time on           
                                    % reshapes and permutes.
    loopTo = k-1;                   % Will only have to loop up to 
                                    % dimension k-1 below.
else                                % If x is a matrix, have to loop over 
    loopTo = k;                     % all dimensions.
end                                 

if k>2 || m>1                       % Now loop over remaining dimensions,  
    x = reshape(x,xsiz);            % inf any. Reshape x into an array of 
    for i=2:loopTo                  % dimension R(1),C(2),...,C(k)or R(k),m.
        if comp(i)                  % If Q{i} is not identity: Create
            dims = 1:k+1;           % vector to re-shuffle dimensions.
            dims(i) = [];           % Put dimension i first (by permute),
            dims = [i, dims];       %#ok<AGROW> % e.g. order [2,1,3,4,5]  
                                    % for i=2 and k=4. Turn off Matlab's
                                    % warning for size change.
            Xmat = reshape( permute(x,dims), [C(i), prod(xsiz)/C(i)] );
                                    % Then bring array into matrix with
            Xmat = Q{i}*Xmat;       % N(i) rows, ex: N(2)-by-N(1)*N(3)*...
                                    % *N(4)*m and multiply by Q{i}.
            xsiz(i) = R(i);         % Changed dimensionality of x.
            x = ipermute( reshape(Xmat,[R(i), xsiz(dims(2:k+1))]), dims );
            %disp(size(Q{i}))
        end                         % Reshape back to array, ex: to dim.
    end                             % N(2),N(1),N(3),N(4),m, and inverse-
                                    % permute to go back to orginal array,
end                                 % ex: dim. N(1),N(2),N(3),N(4),m.
                                    
x = reshape(x,[prod(R),m]);         % Then give back result as matrix.
