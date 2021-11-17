function test_skewHamildeflZ_nb
% Function for testing the routine ZGHUDP (block code). Random matrices of all
% orders m = 1:20 are used and the relative errors are checked.

%
%   Contributor:
%   V. Sima, Dec. 2010.
%
%   Revisions:
%   V. Sima, July 2013, July 2014.
%   M. Voigt, July 2013.
%

disp(' ')
disp('Random tests of the routine ZGHUDP')
disp(' ')

count = 0;
tol   = 100*eps^( 2/3 );
tole  = sqrt( eps );  tolm = 10^5*eps;
max_err = 0.0;  failures = 0;

cmpq0 = 0;  cmpq1 = 1;
orth0 = 0;  orth1 = 1;  orth2 = 2;

% Check default values.

m = 4;  n = 2*m;
nb = 2;
%
A = rand( m ) + rand( m )*1i;  DE = rand( m, m+1 ) + rand( m, m+1 )*1i;  
B = rand( m ) + rand( m )*1i;  FG = rand( m, m+1 ) + rand( m, m+1 )*1i;
%
[ Alpha,  Beta,  Q,  neig  ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1, orth0 );
[ Alpha1, Beta1, Q1, neig1 ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1, orth1 );
[ Alpha2, Beta2, Q2, neig2 ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1, orth2 );
[ Alpha3, Beta3, Q3, neig3 ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1 );
[ Alpha4, Beta4            ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq0 );
[ Alpha5, Beta5            ] = skewHamildeflZ_nb( A, DE, B, FG, nb );
[ Alpha6, Beta6            ] = skewHamildeflZ_nb( A, DE, B, FG );
%
err = max( [ norm( Alpha - Alpha1, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha2, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha3, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha4, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha5, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha6, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Beta  - Beta1,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta2,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta3,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta4,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta5,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta6,  1 )/max( 1, norm( Beta,  1 ) ) ] );
%
err = max( [ neig - neig1, neig - neig2, neig - neig3, err ] );
%
erq = norm( Q - Q3 )/n;
if err ~= 0 || erq > eps, 
    disp( 'Check default values:' )
    if err > tol || erq > tole,  
        disp( 'Failed 1' ),
    else
        disp( 'The most accurate tests (err = 0, erq <= eps) are not satisfied' ) 
        disp( [ 'But err = ', num2str( err), ' and erq = ', num2str( erq) ] )
    end
    disp( ' ' ) 
end
%
% Check functionality. A large m is also used.
%
% Using rand.
%
for m = [ 0 : 20, 210 ]
    nbmax = min( max( m, 1 ), 20 );  
    n = 2*m;
    A = rand( m ) + rand( m )*1i;  DE = rand( m, m+1 ) + rand( m, m+1 )*1i;  
    B = rand( m ) + rand( m )*1i;  FG = rand( m, m+1 ) + rand( m, m+1 )*1i;
    v = diag( DE(:,1:m   ) );  DE(1 : m+1:m*m    ) = ( v - conj( v ) )/2;
    v = diag( DE(:,2:m+1 ) );  DE(m+1:m+1:m*(m+1)) = ( v - conj( v ) )/2;
    v = diag( FG(:,1:m   ) );  FG(1 : m+1:m*m    ) = ( v + conj( v ) )/2;
    v = diag( FG(:,2:m+1 ) );  FG(m+1:m+1:m*(m+1)) = ( v + conj( v ) )/2;
    D = triu( DE(:,2:end) );   D = D - triu( D,  1 )';  
    F = triu( FG(:,2:end) );   F = F + triu( F,  1 )';
    E = tril( DE(:,1:m) );     E = E - tril( E, -1 )';
    G = tril( FG(:,1:m) );     G = G + tril( G, -1 )';
    %
    S = [ A D; E  A' ];
    H = [ B F; G -B' ];
    J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
    %
    if m == 0,  S = [ ];  H = [ ];  end
    evm = eig( H, S );
    %
    % For checking the deflating subspace.
    %
    [ AA, BB, Qq, Z ] = qz( H, S );
    Ev = ordeig( AA, BB );
    for ii = 1 : n,
        if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
            Ev(ii) = imag( Ev(ii) )*1i;
        end
    end
    SELECT = real( Ev < 0 );  ns = numel( find( SELECT == 1 ) );
    [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Z, SELECT );
    %
    for nb = 0 : nbmax,
        for orthm = 0 : 2,
            count = count + 1;
            [ Alpha, Beta, Q, neig ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1, orthm );
            eva = Alpha./Beta;  eva = eva(:);
            err = norm( evm - cmpoles( evm, eva ) )/max( 1, norm( evm ) );
            max_err = max( [ max_err, err ] );
            if max_err > max( 1, n )*tol,  disp( 'Failed 2' ),  return,  end
            %
            Vx = [ Q Zs(:,1:ns) ];  tolr = tolm*n*max( norm( Vx ), 1 );
            err = max( [ abs( neig - ns ) abs( ns - rank( Vx, tolr ) ) ] );
            if err > 0,  
                disp( [ 'Failed 3, orthm = ', num2str( orthm ) ] )
                failures = failures + 1;
            end
        end
    end
end

if failures ~= 0,  disp( ' ' ),  end
if max_err < tole,
    disp( [ 'ZGHUDP :    passed  -- maximum relative error max_err = ', num2str( max_err  ) ] )
    disp( [ '            Number of problems solved using rand      = ', num2str( count    ) ] )
    disp( ' ' )
else
    disp( [ 'ZGHUDP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
    disp( [ '            Number of problems solved using rand      = ', num2str( count   ) ] )
end

if failures > 0,  
    disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
    disp( ' ' )
end        
%
count   = 0;
max_err = 0.0;
%
% Using randn.
%
for m = [ 0 : 20, 210 ],
    nbmax = min( max( m, 1 ), 20 );  
    n = 2*m;
    A = randn( m ) + randn( m )*1i;  DE = randn( m, m+1 ) + randn( m, m+1 )*1i;  
    B = randn( m ) + randn( m )*1i;  FG = randn( m, m+1 ) + randn( m, m+1 )*1i;
    v = diag( DE(:,1:m   ) );  DE(1 : m+1:m*m    ) = ( v - conj( v ) )/2;
    v = diag( DE(:,2:m+1 ) );  DE(m+1:m+1:m*(m+1)) = ( v - conj( v ) )/2;
    v = diag( FG(:,1:m   ) );  FG(1 : m+1:m*m    ) = ( v + conj( v ) )/2;
    v = diag( FG(:,2:m+1 ) );  FG(m+1:m+1:m*(m+1)) = ( v + conj( v ) )/2;
    D = triu( DE(:,2:end) );   D = D - triu( D,  1 )';  
    F = triu( FG(:,2:end) );   F = F + triu( F,  1 )';
    E = tril( DE(:,1:m) );     E = E - tril( E, -1 )';
    G = tril( FG(:,1:m) );     G = G + tril( G, -1 )';
    %
    S = [ A D; E  A' ];
    H = [ B F; G -B' ];
    J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
    %
    if m == 0,  S = [ ];  H = [ ];  end
    evm = eig( H, S );
    %
    % For checking the deflating subspace.
    %
    [ AA, BB, Qq, Z ] = qz( H, S );
    Ev = ordeig( AA, BB );
    for ii = 1 : n,
        if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
            Ev(ii) = imag( Ev(ii) )*1i;
        end
    end
    SELECT = real( Ev < 0 );  ns = numel( find( SELECT == 1 ) );
    [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Z, SELECT );
    %
    for nb = 0 : nbmax,
        for orthm = 0 : 2,
            count = count + 1;
            [ Alpha, Beta, Q, neig ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1, orthm );
            eva = Alpha./Beta;  eva = eva(:);
            err = norm( evm - cmpoles( evm, eva ) )/max( 1, norm( evm ) );
            max_err = max( [ max_err, err ] );
            if max_err > max( 1, n )*tol,  disp( 'Failed 4' ),  return,  end
            %
            Vx = [ Q Zs(:,1:ns) ];  tolr = tolm*n*max( norm( Vx ), 1 );
            err = max( [ abs( neig - ns ) abs( ns - rank( Vx, tolr ) ) ] );
            if err > 0,  
                disp( [ 'Failed 5, orthm = ', num2str( orthm ) ] )
                failures = failures + 1;
            end
        end
    end
end

if failures ~= 0,  disp( ' ' ),  end
if max_err < tole,
    disp( [ 'ZGHUDP :    passed  -- maximum relative error max_err = ', num2str( max_err  ) ] )
    disp( [ '            Number of problems solved using randn     = ', num2str( count    ) ] )
    disp( ' ' )
else
    disp( [ 'ZGHUDP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
    disp( [ '            Number of problems solved using randn     = ', num2str( count   ) ] )
end

if failures > 0,  
    disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
    disp( ' ' )
end        

