function perf_skewHamildeflZ_nb
% Function for testing the performance of the routine ZGHUDP (block code).
% Random matrices of orders 2*list_m are used and the relative errors are checked.
% For each m in list_m, nrtest are performed, and the mean values are
% returned.  Default values: list_m = 10 : 5 : 20;  nrtest = 5;

%
%   Contributor:
%   V. Sima, Dec. 2010.
%
%   Revisions:
%   V. Sima, Jul. 2013.
%   M. Voigt, Jul. 2013.
%

disp(' ')
disp( 'Random performance tests of the routine ZGHUDP' )
disp(' ')

if ~exist( 'list_nb', 'var' ) || isempty( list_nb ),  list_nb = 32;  end
nrnb = numel( list_nb );

if ~exist( 'list_m', 'var' ) || isempty( list_m ),  list_m = 10 : 5 : 20;  end
if ~exist( 'nrtest', 'var' ) || isempty( nrtest ),  nrtest = 5;  end
nlistm = numel( list_m );

tic

counter = 0;
tol     = 100*eps^( 2/3 );
tole    = sqrt( eps );  tolm = 10^5*eps;
max_err = 0.0;  failures = 0;

cmpq1 = 1;  orthm = 2;  % Only SVD is used.

err      = zeros( nrtest, nrnb );
ert      = zeros( nrtest, nrnb );
timing   = zeros( nrtest, nrnb+1 );
mtiming  = zeros( nlistm, nrnb+1 );
meanerr  = zeros( nlistm, nrnb );
meanert  = zeros( nlistm, nrnb );
%
% Using rand.
%
for m = list_m,
    n = 2*m;
    counter = counter + 1;
    disp( [ 'n = ', num2str( n ), ',  m = ', num2str( m ) ] )
    for j = 1 : nrtest,
        %
        % Generate the input matrices.
        %
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
        clear D E F G v
        %
        if m == 0,  S = [ ];  H = [ ];  end
        evm = eig( H, S );
        %
        % MATLAB computation.
        %
        time = tic;
        [ AA, BB, Qq, Z ] = qz( H, S );
        Ev = ordeig( AA, BB );
        for jj = 1 : n,
            if abs( real( Ev(jj) ) ) / abs( Ev(jj) ) < tolm*n,
                Ev(jj) = imag( Ev(jj) )*1i;
            end
        end
        SELECT = real( Ev < 0 );  ns = numel( find( SELECT == 1 ) );
        [ ~, ~, ~, Zs ] = ordqz( AA, BB, Qq, Z, SELECT );
        time = toc( time );
        timing(j,1) = time;
        clear AA BB Ev H Qq S SELECT Z
        %
        % FORTRAN computation.
        %
        k = 0;
        for nb = list_nb,
            k = k + 1;
            time = tic;
            [ Alpha, Beta, Q, neig ] = skewHamildeflZ_nb( A, DE, B, FG, nb, cmpq1, orthm );
            ev = Alpha./Beta;  ev = ev(:);
            time = toc( time );
            timing(j,1+k) = time;
            err(j,k) = norm( evm - cmpoles( evm, ev ) )/max( 1, norm( evm ) );
            max_err = max( [ max_err, err(j,k) ] );
            if max_err > max( 1, n )*tol,  disp( 'Failed 1' ),  return,  end
            %
            Vx = [ Q Zs(:,1:ns) ];  tolr = tolm*n*max( norm( Vx ), 1 );
            ert(j,k) = max( [ abs( neig - ns ) abs( ns - rank( Vx, tolr ) ) ] );
            if ert(j,k) > 0,  
                if failures < 2,  disp( 'Failed 2' ),  end
                failures = failures + 1;
            end
        end
    end
    if nrtest > 1,
        mtiming (counter,:) = mean( timing  );
        meanerr (counter,:) = mean( err     );
        meanert (counter,:) = mean( ert     );
    else
        mtiming (counter,:) = timing;
        meanerr (counter,:) = err;
        meanert (counter,:) = ert;
    end
    %
    disp( ' ' )
    disp( 'Mean relative errors for eigenvalues and deflating subspace computation with ZGHUDP' )
    disp( ' ' )
    disp( '------------------------------------------' )
    disp( '   nb    eigenvalues   deflating subspace ' )
    disp( '------------------------------------------' )
    for jj = 1 : nrnb,
        disp( [ sprintf( '%5d', list_nb(jj) ), ...
                sprintf( '%15.5g', meanerr(counter,jj) ), ...
                sprintf( '%15.5g', meanert(counter,jj) ) ] )
    end
    disp( '------------------------------------------' )
    disp( ' ' )
    %
    disp( 'Timings (sec) for deflating subspace computation' )
    disp( ' ' )
    disp( '-------------------' )
    disp( '        eig+ordeig ' )
    disp( '-------------------' )
    disp( [ '     ', sprintf( '%13.5g', mtiming (counter,1) ) ] )
    disp( '-------------------' )
    disp( '   nb      ZGHUDP  ' )
    disp( '-------------------' )
    for jj = 1 : nrnb,
        disp( [ sprintf( '%5d', list_nb(jj) ), ...
                sprintf( '%13.5g', mtiming (counter,jj+1) ) ] )
    end
    disp( '-------------------' )
    disp( ' ' )
end

toc

disp( [ 'Sum of timings = ', sprintf( '%13.5g', nrtest*sum( sum( mtiming ) ) ) ] )
disp( ' ' )

countT = counter*nrtest;
if max_err < tole,
   disp( [ 'ZGHUDP :    passed  -- maximum relative error max_err = ', num2str( max_err  ) ] )
   disp( [ '            Number of problems solved using rand      = ', num2str( countT    ) ] )
   disp( ' ' )
else
   disp( [ 'ZGHUDP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using rand      = ', num2str( countT   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
   disp( ' ' )
end        

currdate = date;
filnam = [ 'perf_skewHamildeflZ_' currdate '_large' ];
%
existfile = exist( [ filnam '.mat' ] ) == 2;
while existfile,
    answ = input( [ 'File ', [ filnam '.mat' ], ' already exists. Overwrite it (y/n)? ' ], 's' );
    if strcmp( lower( answ(1) ), 'n' ),
        filnam = input( [ 'New file name (without mat extension): ' ], 's' );
        existfile = exist( [ filnam '.mat' ] ) == 2;
    else
        existfile = 0;
    end
end
%
save( filnam, 'counter', 'list_m', 'list_nb', 'meanerr', 'meanert', 'mtiming', 'nrtest' )
%
disp( [ 'Main results saved in the file ', filnam ] )
