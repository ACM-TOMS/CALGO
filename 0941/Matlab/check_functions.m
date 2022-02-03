function check_functions()
%CHECK_FUNCTIONS Check of main toolbox functionality
%
% CHECK_FUNCTIONS() checks the output of the main functions of the toolbox
% for random tensors of varying size and varying ranks.
% Warning: This may take a long time.
%
% See also TEST_FUNCTIONS, TEST_ERRORS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Tolerance to account for roundoff error
tol = 10^(-13);

randn('state',0);
rand('state',0);

% Tensors of order 2,3,4,5 are tested.
% Specify the sizes to be tested.
sizes{2} = [1,1;
            10,1;
            1,10;
            5,6;
            7,7];
sizes{3} = [1,1,1;
            10,1,1;
            1,1,10;
            5,6,4;
            7,7,7];
sizes{4} = [1,1,1,1;
            10,1,1,1;
            1,1,1,10;
            5,6,4,7;
            7,7,7,7];
sizes{5} = [1,1,1,1,1;
            10,1,1,1,1;
            1,1,1,1,10;
            5,6,4,7,3;
            7,7,7,7,7];
treetypes = {'first_separate','first_pair_separate','TT','balanced'};


% ---------------------------------------------
% Tests of exact operations with random tensors
%----------------------------------------------

tests = 0;
suctests = 0;
for d = 2:5,
    szs = sizes{d};
    for szind = 1:size(szs,1)
        sz = szs( szind,: );
        for treeind = 1:length(treetypes),
            treetype = treetypes{treeind};
            for rankXtype = 1:3,
                if rankXtype == 1,
                    x = htenrandn(sz, '', ones(1,2*d-1), treetype );
                elseif rankXtype == 2,
                    x = htenrandn(sz, '', 3*ones(1,2*d-1), treetype );
                else
                    x = htenrandn(sz, '', [], treetype );
                end
                normX = max(1,norm(x));
                
                % ---------------------------------------------------
                % Test functions that involve a single htucker tensor
                %----------------------------------------------------

                for ind = 1:2*d-1,
                   tests = tests + 1;
                   if check_change_root(x,ind,tol*normX)
                      suctests = suctests + 1;
                   else
                      warning('Test %i failed: change_root.',tests);
                   end
                end
                
                tests = tests + 1;
                if check_gramians(orthog(x),tol*normX^2)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: gramians.',tests);
                end
                
                tests = tests + 1;
                if check_ipermute(x,randperm(d),tol*normX)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: ipermute.',tests);
                end

                tests = tests + 1;
                if check_mrdivide(x,pi,tol*normX/pi)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: mrdivide.',tests);
                end
                
                tests = tests + 1;
                if check_mtimes(x,pi,tol*normX*pi)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: mtimes.',tests);
                end

                tests = tests + 1;
                if check_norm(x,tol*normX)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: norm.',tests);
                end
                
                if sz(end)~=1,
                    tests = tests + 1;
                    y = randn(sz);  normY = max(1,norm(y(:)));
                    if check_norm_diff(x,y,tol*max(normX,normY))
                        suctests = suctests + 1;
                    else
                        warning('Test %i failed: norm_diff.',tests);
                    end
                end
                
                tests = tests + 1;
                if check_orthog(x,tol*normX)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: orthog.',tests);
                end

                tests = tests + 1;
                if check_permute(x,randperm(d),tol*normX)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: permute.',tests);
                end
                
                tests = tests + 1;
                if check_power(x,tol*normX^2)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: power.',tests);
                end
                
                tests = tests + 1;
                if check_singular_values(x,sqrt(tol)*normX)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: singular_values.',tests);
                    return
                end

                A = [];
                normA = 1;
                for j = 1:d,
                    A{j} = rand(5,sz(j));
                    normA = max(normA,norm(A{j}));
                end
                tests = tests + 1;
                if check_ttm(x,A,tol*normX*normA)
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: ttm.',tests);
                end
                
                for j = 1:d,
                    tests = tests + 1;
                    if check_ttmi(x,A{j},j,tol*normX*normA)
                        suctests = suctests + 1;
                    else
                        warning('Test %i failed: ttm, individual.',tests);
                    end
                end

                A = [];
                normA = 1;
                for j = 1:d,
                    A{j} = rand(sz(j),5);
                    normA = max(normA,norm(A{j}));
                end
                tests = tests + 1;
                if check_ttm(x,A,tol*normX*normA,'t')
                    suctests = suctests + 1;
                else
                    warning('Test %i failed: ttm, transpose.',tests);
                end
                
                for j = 1:d,
                    tests = tests + 1;
                    if check_ttmi(x,A{j},j,tol*normX*normA,'t')
                        suctests = suctests + 1;
                    else
                        warning('Test %i failed: ttm, individual, transpose.',tests);
                    end
                end
                
                tests = tests + 1;
                v = [];
                for j = 1:d,
                    v{j} = randn(sz(j),1);
                    v{j} = v{j} / norm( v{j} );
                end
                if check_ttv(x,v,tol*normX)
                   suctests = suctests + 1;
                else
                   warning('Test %i failed: ttv.',tests);
                end

                
                % -----------------------------------------------
                % Test functions that involve two htucker tensors
                % -----------------------------------------------
                for rankYtype = 1:3,
                    if rankYtype == 1,
                        y = htenrandn(sz, '', ones(1,2*d-1), treetype );
                    elseif rankYtype == 2,
                        y = htenrandn(sz, '', 3*ones(1,2*d-1), treetype );
                    else
                        y = htenrandn(sz, '', [], treetype );
                    end
                    normY = max(1,norm(y));

                    tests = tests + 1;
                    if check_innerprod(x,y,tol*normX*normY)
                       suctests = suctests + 1;
                    else
                        warning('Test %i failed: innerprod.',tests);
                    end
                    
                    tests = tests + 1;
                    if check_minus(x,y,tol*max(normX,normY))
                       suctests = suctests + 1;
                    else
                        warning('Test %i failed: minus.',tests);
                    end
                    
                    tests = tests + 1;
                    if check_plus(x,y,tol*max(normX,normY))
                       suctests = suctests + 1;
                    else
                        warning('Test %i failed: plus.',tests);
                    end
                    
                    tests = tests + 1;
                    if check_times(x,y,tol*normX*normY)
                       suctests = suctests + 1;
                    else
                        warning('Test %i failed: times.',tests);
                    end
                    
                    tests = tests + 1;
                    dims = 1:ceil(d/2);
                    if check_ttt(x,y,dims,tol*normX*normY)
                       suctests = suctests + 1;
                    else
                        warning('Test %i failed: ttt.',tests);
                    end
                    
                end
            end
            fprintf('%i of %i tests successfully passed.\n',suctests, tests );
        end
    end
end

% ----------------------------------------------------------
% Tests of truncating operations with function-based tensors
%-----------------------------------------------------------

for d = 2:6,
   opts.max_rank = 50; opts.abs_eps = 10^(-5);
   A = reciproc_sum(d, 100, 0.01, 1, 50); 
   x = htensor.truncate_cp(A,opts);
   for tol = [10^(-4), 10^(-3), 10^(-2), 10^(-1)],
       opts.abs_eps = tol;
       
       tests = tests + 1;
       y = truncate_std(x, opts);
       if norm(orthog(x-y)) <= tol,
           suctests = suctests + 1;
       else
          warning('Test %i failed: truncate_std.',tests); 
       end

       tests = tests + 1;
       y = truncate_nonorthog(x, opts);
       if norm(orthog(x-y)) <= tol,
           suctests = suctests + 1;
       else
          warning('Test %i failed: truncate_nonorthog.',tests); 
       end
   end
   fprintf('%i of %i tests successfully passed.\n',suctests, tests );
end

for d = 2:4,
   opts.max_rank = 50;
   if d==4,
       n = 10;
   else
       n = 50;
   end
   A = reciproc_sum(d, n, 0.01, 1); 
   for tol = [10^(-4), 10^(-3), 10^(-2), 10^(-1)],
       opts.abs_eps = tol;
       
       tests = tests + 1;
       x = htensor.truncate_ltr(A, opts);
       if norm( x(:) - A(:) ) <= tol,
           suctests = suctests + 1;
       else
          warning('Test %i failed: truncate_ltr.',tests); 
       end

       tests = tests + 1;
       x = htensor.truncate_rtl(A, opts);
       if norm( x(:) - A(:) ) <= tol,
           suctests = suctests + 1;
       else
          warning('Test %i failed: truncate_rtl.',tests); 
       end
   end
   fprintf('%i of %i tests successfully passed.\n',suctests, tests );
end

function res = check_change_root(x,ind,tol)
% Check that changing the root does not alter the tensor
y = change_root(x,ind);
err = norm( x(:) - y(:) );
res = ( err <= tol );

function res = check_gramians(x,tol)
% Check Gramians by comparing their singular values with the matricization
err = 0;
G = gramians(x);
for j = 2:length(G),
   A = matricize(full(x),x.dims{j});
   sA = svd(A).^2;
   sG = svd(G{j});
   l = min(length(sA),length(sG));
   err = max( err, norm( sG(1:l) - sA(1:l) ) );
end
res = ( err <= tol );

function res = check_innerprod(x,y,tol)
% Compare htucker inner product with inner product of vectors
err = abs( innerprod(x,y) - dot( x(:), y(:) ) );
res = ( err <= tol );

function res = check_ipermute(x,perm,tol)
% Compare inverse permutation of htucker tensor with inverse permutation of array
y = ipermute(x,perm);
z = ipermute(full(x),perm);
err = norm( z(:) - y(:) );
res = ( err <= tol );

function res = check_minus(x,y,tol)
% Compare htucker subtraction with subtraction of vectors
z = minus(x,y);
err = norm( z(:) - ( x(:) - y(:) ) );
res = ( err <= tol );

function res = check_mrdivide(x,scalar,tol)
% Compare htucker scalar division with scalar division of vector
y = x / scalar;
err = norm( x(:) / scalar - y(:) );
res = ( err <= tol );

function res = check_mtimes(x,scalar,tol)
% Compare htucker scalar multiplication with scalar multiplication of vector
y = scalar*x;
err = norm( scalar*x(:) - y(:) );
res = ( err <= tol );

function res = check_norm(x,tol)
% Compare htucker norm with norm of vector
err = abs( norm(x) - norm( x(:) ) );
res = ( err <= tol );

function res = check_norm_diff(x,y,tol)
% Compare htucker norm difference with norm difference of vectors
err = abs( norm_diff(x,y) - norm( x(:) - y(:) ) );
res = ( err <= tol );

function res = check_orthog(x,tol)
% Compare tensor before and after orthogonalization
y = orthog(x);
err = norm( x(:) - y(:) );
res = ( err <= tol );

function res = check_power(x,tol)
% Compare htucker elementwise square with elementwise square of vector
y = power(x,2);
err = norm( x(:).^2 - y(:) );
res = ( err <= tol );

function res = check_permute(x,perm,tol)
% Compare permutation of htucker tensor with permutation of array
y = permute(x,perm);
z = permute(full(x),perm);
err = norm( z(:) - y(:) );
res = ( err <= tol );

function res = check_plus(x,y,tol)
% Compare htucker addition with addition of vectors
z = plus(x,y);
err = norm( z(:) - ( x(:) + y(:) ) );
res = ( err <= tol );

function res = check_singular_values(x,tol)
% Check singular values by comparing with the matricization
err = 0;
s = singular_values(x);
for j = 2:length(s),
   A = matricize(full(x),x.dims{j});
   sA = svd(A);
   sG = s{j};
   l = min(length(sA),length(sG));
   err = max( err, norm( sA(1:l) - sG(1:l) ) );
end
res = ( err <= tol );

function res = check_times(x,y,tol)
% Compare htucker elementwise product with elementwise product of vectors
z = times(x,y);
err = norm( z(:) - ( x(:) .* y(:) ) );
res = ( err <= tol );

function res = check_ttt(x,y,dims,tol)
% Compare htucker tensor-times-tensor with full tensor-times-tensor
z  = ttt(x,y,dims);
zf = ttt(full(x),full(y),dims);
err = norm( z(:) - zf(:) );
res = ( err <= tol );

function res = check_ttm(x,A,tol,transp)
% Compare htucker N-mode matrix products with array N-mode matrix products
if nargin < 4,
    y = ttm(x,A);
    z = ttm(full(x),A);
else
    y = ttm(x,A,transp);
    z = ttm(full(x),A,transp);    
end
err = norm( z(:) - y(:) );
res = ( err <= tol );

function res = check_ttmi(x,A,j,tol,transp)
% Compare htucker N-mode matrix products with array N-mode matrix products
if nargin < 5,
    y = ttm(x,A,j);
    z = ttm(full(x),A,j);
else
    y = ttm(x,A,j,transp);
    z = ttm(full(x),A,j,transp);
end
err = norm( z(:) - y(:) );
res = ( err <= tol );

function res = check_ttv(x,v,tol)
% Compare htucker ttv with correpsonding array N-mode matrix products
y = ttv(x,v);
z = ttm(full(x),v,'t');;
err = abs( z(:) - y(:) );
res = ( err <= tol );
