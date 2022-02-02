function colamd_test (spumoni)
%
% colamd_test (spumoni)
%
% COLAMD and SYMAMD testing function.  Here we try to give colamd and symamd
% every possible type of matrix and erroneous input that they may encounter. 
% We want either a valid permutation returned or we want them to fail
% gracefully.  If the optional spumoni argument is present, additional
% diagnostic messages are printed during the test:
%
% spumoni:
%   0: very little, just a progress meter (default)
%   1: lots
%   2: far too much
%
% You are prompted as to whether or not the colamd and symand routines and
% the test mexFunctions are to be compiled.
%
% Tim Davis
% September 8, 2003.  Version 2.3.
% http://www.cise.ufl.edu/research/sparse/colamd/

help colamd_test

s = input (...
'Compile colamd, symand, and the test codes? (y/n, default is yes): ', 's') ;

do_compile = 1 ;
if (~isempty (s))
    if (s (1) == 'n' | s (1) == 'N')
	do_compile = 0 ;
    end
end

if (do_compile)
    fprintf ('Compiling colamd, symamd, and test mexFunctions.\n') ;
    fprintf ('mex -O colamdmex.c colamd.c\n') ;
    mex -O colamdmex.c colamd.c
    fprintf ('mex -O symamdmex.c colamd.c\n') ;
    mex -O symamdmex.c colamd.c
    fprintf ('mex -O colamdtestmex.c colamd.c\n') ;
    mex -O colamdtestmex.c colamd.c
    fprintf ('mex -O symamdtestmex.c colamd.c\n') ;
    mex -O symamdtestmex.c colamd.c
    fprintf ('Done compiling.\n') ; 
end
    
fprintf ('\nThe following codes will be tested:\n') ;
which colamd 
which symamd
which colamdmex
which symamdmex

fprintf ('\nStarting the tests.  Please be patient.\n') ;

rand ('state', 0) ;
randn ('state', 0) ;
if (nargin < 1)
    spumoni = 0 ;
end
old = spparms ('spumoni') ;
spparms ('spumoni', spumoni) ;


fprintf ('Null matrices') ;
A = zeros (0,0) ;
A = sparse (A) ;

[p, stats] = colamd (A, [.5 .5 spumoni]) ;
check_perm (p, A) ;

[p, stats] = symamd (A, [.5 spumoni]) ;
check_perm (p, A) ;

A = zeros (0, 100) ;
A = sparse (A) ;
[p, stats] = colamd (A, [.5 .5 spumoni]) ;
check_perm (p, A) ;

A = zeros (100, 0) ;
A = sparse (A) ;
[p, stats] = colamd (A, [.5 .5 spumoni]) ;
check_perm (p, A) ;
fprintf (' OK\n') ;


fprintf ('Matrices with a few dense row/cols\n') ;

for trial = 1:20

    % random square unsymmetric matrix
    A = rand_matrix (1000, 1000, 1, 10, 20) ;
    [m n] = size (A) ;


    for tol = [0 0.001 0.005 0.01:0.01:0.10 0.5:.1:1]

	[p, stats] = colamd (A, [tol tol spumoni]) ;
	check_perm (p, A) ;

	B = A + A' ;
	[p, stats] = symamd (B, [tol spumoni]) ;
	check_perm (p, A) ;

	[p, stats] = colamd (A, [tol 1 spumoni]) ;
	check_perm (p, A) ;

	[p, stats] = colamd (A, [1 tol spumoni]) ;
	check_perm (p, A) ;

	fprintf ('.') ;

    end
end
fprintf (' OK\n') ;

fprintf ('General matrices\n') ;
for trial = 1:400

    % matrix of random mtype
    mtype = irand (3) ;
    A = rand_matrix (2000, 2000, mtype, 0, 0) ;
    [m n] = size (A) ;
    p = colamd (A) ;
    check_perm (p, A) ;
    if (mtype == 3)
	p = symamd (A) ;
	check_perm (p, A) ;
    end

    fprintf ('.') ;
end
fprintf (' OK\n') ;

fprintf ('Test error handling with invalid inputs\n') ;

% Check different erroneous input.
for trial = 1:30

    A = rand_matrix (1000, 1000, 2, 0, 0) ;
    [m n] = size (A) ;

    for err = 1:13

        p = Tcolamd (A, [1 1 spumoni 0 err]) ;
        if p ~= -1 
	    check_perm (p, A) ;
	end

	if (err == 1)
	    % check different (valid) input args to colamd
	    p = Acolamd (A, spumoni) ;

	    p2 = Acolamd (A, spumoni, [.5 .5 spumoni 0 0]) ;
	    if (any (p ~= p2))
		error ('colamd: mismatch 1!') ;
	    end
	    [p2 stats] = Acolamd (A, spumoni) ;
	    if (any (p ~= p2))
		error ('colamd: mismatch 2!') ;
	    end
	    [p2 stats] = Acolamd (A, spumoni, [.5 .5 spumoni 0 0]) ;
	    if (any (p ~= p2))
		error ('colamd: mismatch 3!') ;
	    end
	end

	B = A'*A ;
        p = Tsymamd (B, [1 spumoni err]) ;
        if p ~= -1 
	    check_perm (p, A) ;
	end

	if (err == 1)

	    % check different (valid) input args to symamd
	    p = Asymamd (B, spumoni) ;
	    check_perm (p, A) ;
	    p2 = Asymamd (B, spumoni, [.5 spumoni 0]) ;
	    if (any (p ~= p2))
		error ('symamd: mismatch 1!') ;
	    end
	    [p2 stats] = Asymamd (B, spumoni) ;
	    if (any (p ~= p2))
		error ('symamd: mismatch 2!') ;
	    end
	    [p2 stats] = Asymamd (B, spumoni, [.5 spumoni 0]) ;
	    if (any (p ~= p2))
		error ('symamd: mismatch 3!') ;
	    end
	end

	fprintf ('.') ;
    end

end
fprintf (' OK\n') ;

fprintf ('Matrices with a few empty columns\n') ;

for trial = 1:400

    % some are square, some are rectangular
    n = 0 ;
    while (n < 5)
	A = rand_matrix (1000, 1000, irand (2), 0, 0) ;
	[m n] = size (A) ;
    end

    % Add 5 null columns at random locations.
    null_col = randperm (n) ;
    null_col = sort (null_col (1:5)) ;
    A (:, null_col) = 0 ;

    % Order the matrix and make sure that the null columns are ordered last.
    [p, stats] = colamd (A, [1 1 spumoni]) ;
    check_perm (p, A) ;

    if (stats (2) ~= 5)
	error ('colamd: wrong number of null columns') ;
    end
    if (any (null_col ~= p ((n-4):n)))
	error ('colamd: Null cols are not ordered last in natural order') ;
    end

    fprintf ('.') ;

end
fprintf (' OK\n') ;

fprintf ('Matrices with a few empty rows and columns\n') ;

for trial = 1:400

    % symmetric matrices
    n = 0 ;
    while (n < 5)
	A = rand_matrix (1000, 1000, 3, 0, 0) ;
	[m n] = size (A) ;
    end

    % Add 5 null columns and rows at random locations.
    null_col = randperm (n) ;
    null_col = sort (null_col (1:5)) ;
    A (:, null_col) = 0 ;
    A (null_col, :) = 0 ;

    % Order the matrix and make sure that the null rows/cols are ordered last.
    [p,stats] = symamd (A, [1 spumoni]) ;
    check_perm (p, A) ;

    % find actual number of null rows and columns
    Alo = tril (A, -1) ;
    nnull = length (find (sum (Alo') == 0 & sum (Alo) == 0)) ;

    if (stats (2) ~= nnull | nnull < 5)
	error ('symamd: wrong number of null columns') ;
    end
    if (any (null_col ~= p ((n-4):n)))
	error ('symamd: Null cols are not ordered last in natural order') ;
    end

    fprintf ('.') ;

end
fprintf (' OK\n') ;

fprintf ('Matrices with a few empty rows\n') ;

% Test matrices with null rows inserted.

for trial = 1:400

    m = 0 ;
    while (m < 5)
	A = rand_matrix (1000, 1000, 2, 0, 0) ;
	[m n] = size (A) ;
    end

    % Add 5 null rows at random locations.
    null_row = randperm (m) ;
    null_row = sort (null_row (1:5)) ;
    A (null_row, :) = 0 ;

    p = colamd (A, [.5 .5 spumoni]) ;
    check_perm (p, A) ;
    if (stats (1) ~= 5)
	error ('colamd: wrong number of null rows') ;
    end
    fprintf ('.') ;
end
fprintf (' OK\n') ;


fprintf ('\ncolamd and symamd:  all tests passed\n\n') ;
spparms ('spumoni', old) ;

%-------------------------------------------------------------------------------
% Acolamd:  compare colamd and Tcolamd results
%-------------------------------------------------------------------------------

function [p,stats] = Acolamd (S, spumoni, knobs)

if (nargin < 2)
    spumoni = 0 ;
end

if (nargin < 3)
    if (nargout == 1)
	[p] = colamd (S) ;
	[p1] = Tcolamd (S, [.5 .5 spumoni 0 0]) ;
    else
	[p, stats] = colamd (S) ;
	[p1, stats1] = Tcolamd (S, [.5 .5 spumoni 0 0]) ;
    end
else
    if (nargout == 1)
	[p] = colamd (S, knobs (1:3)) ;
	[p1] = Tcolamd (S, knobs) ;
    else
	[p, stats] = colamd (S, knobs (1:3)) ;
	[p1, stats1] = Tcolamd (S, knobs) ;
    end
end

check_perm (p, S) ;
check_perm (p1, S) ;

if (any (p1 ~= p))
    error ('Acolamd mismatch!') ;
end


%-------------------------------------------------------------------------------
% Asymamd:  compare symamd and Tsymamd results
%-------------------------------------------------------------------------------

function [p,stats] = Asymamd (S, spumoni, knobs)

if (nargin < 2)
    spumoni = 0 ;
end

if (nargin < 3)
    if (nargout == 1)
	[p] = symamd (S) ;
	[p1] = Tsymamd (S, [.5 spumoni 0]) ;
    else
	[p, stats] = symamd (S) ;
	[p1, stats1] = Tsymamd (S, [.5 spumoni 0]) ;
    end
else
    if (nargout == 1)
	[p] = symamd (S, knobs (1:2)) ;
	[p1] = Tsymamd (S, knobs) ;
    else
	[p, stats] = symamd (S, knobs (1:2)) ;
	[p1, stats1] = Tsymamd (S, knobs) ;
    end
end

if (any (p1 ~= p))
    error ('Asymamd mismatch!') ;
end


%-------------------------------------------------------------------------------
% check_perm:  check for a valid permutation vector
%-------------------------------------------------------------------------------

function check_perm (p, A)

if (isempty (A) & isempty (p))
    % empty permutation vectors of empty matrices are OK
    return
end

if (isempty (p))
    error ('bad permutation: cannot be empty') ;
end

[m n] = size (A) ;
[pm pn] = size (p) ;
if (pn == 1)
    % force p to be a row vector
    p = p' ;
    [pm pn] = size (p) ;
end

if (n ~= pn)
    error ('bad permutation: wrong size') ;
end

if (pm ~= 1) ;
    % p must be a vector
    error ('bad permutation: not a vector') ;
else
    if (any (sort (p) - (1:pn)))
	error ('bad permutation') ;
    end
end

%-------------------------------------------------------------------------------
% irand: return a random integer between 1 and n
%-------------------------------------------------------------------------------

function i = irand (n)
i = min (n, 1 + floor (rand * n)) ;

%-------------------------------------------------------------------------------
% rand_matrix:  return a random sparse matrix
%-------------------------------------------------------------------------------

function A = rand_matrix (nmax, mmax, mtype, drows, dcols)
%
% A = rand_matrix (nmax, mmax, mtype, drows, dcols)
%
% A binary matrix of random size, at most nmax-by-mmax, with drows dense rows
% and dcols dense columns.
%
% mtype 1: square unsymmetric (mmax is ignored)
% mtype 2: rectangular
% mtype 3: symmetric (mmax is ignored)

n = irand (nmax) ;
if (mtype ~= 2)
    % square
    m = n ;
else
    m = irand (mmax) ;
end

A = sprand (m, n, 10 / max (m,n)) ;

if (drows > 0)
    % add dense rows
    for k = 1:drows
	i = irand (m) ;
	nz = irand (n) ;
	p = randperm (n) ;
	p = p (1:nz) ;
	A (i,p) = 1 ;
    end
end

if (dcols > 0)
    % add dense cols
    for k = 1:dcols
	j = irand (n) ;
	nz = irand (m) ;
	p = randperm (m) ;
	p = p (1:nz) ;
	A (p,j) = 1 ;
    end
end

A = spones (A) ;

% ensure that there are no empty columns
d = find (full (sum (A)) == 0) ;
A (m,d) = 1 ;

% ensure that there are no empty rows
d = find (full (sum (A,2)) == 0) ;
A (d,n) = 1 ;

if (mtype == 3)
    % symmetric
    A = A + A' + speye (n) ;
end

A = spones (A) ;

%-------------------------------------------------------------------------------
% Tcolamd:  run colamd in a testing mode
%-------------------------------------------------------------------------------

function [p,stats] = Tcolamd (S, knobs)

if (nargout <= 1 & nargin == 1)
    p = colamdtestmex (S) ;
elseif (nargout <= 1 & nargin == 2)
    p = colamdtestmex (S, knobs) ;
elseif (nargout == 2 & nargin == 1)
    [p, stats] = colamdtestmex (S) ;
elseif (nargout == 2 & nargin == 2)
    [p, stats] = colamdtestmex (S, knobs) ;
else
    error ('colamd:  incorrect number of input and/or output arguments') ;
end

if (p (1) ~= -1)
    [ignore, q] = sparsfun ('coletree', S (:,p)) ;
    p = p (q) ;
    check_perm (p, S) ;
end

%-------------------------------------------------------------------------------
% Tsymamd: run symamd in a testing mode
%-------------------------------------------------------------------------------

function [p, stats] = Tsymamd (S, knobs)

if (nargout <= 1 & nargin == 1)
    p = symamdtestmex (S) ;
elseif (nargout <= 1 & nargin == 2)
    p = symamdtestmex (S, knobs) ;
elseif (nargout == 2 & nargin == 1)
    [p, stats] = symamdtestmex (S) ;
elseif (nargout == 2 & nargin == 2)
    [p, stats] = symamdtestmex (S, knobs) ;
else
    error ('symamd:  incorrect number of input and/or output arguments') ;
end

if (p (1) ~= -1)
    [ignore, q] = sparsfun ('symetree', S (p,p)) ;
    p = p (q) ;
    check_perm (p, S) ;
end

