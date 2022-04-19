function perf_skewHamildefl_nb
% Function for testing the performance of the routine DGHUDP (block code).
% Random matrices of orders 2*list_m are used and the relative errors are checked.
% For each m in list_m, nrtest are performed, and the mean values are
% returned.  Default values: list_m = 10 : 5 : 20;  nrtest = 5;

%
%   Contributor:
%   V. Sima, Mar. 2011.
%
%   Revisions:
%   V. Sima, Jul. 2013.
%   M. Voigt, Jul. 2013.
%

disp( ' ' )
disp( 'Random performance tests of the routine DGHUDP' )
disp( ' ' )

if ~exist( 'list_nb', 'var' ) || isempty( list_nb ),  list_nb = 32;  end
nrnb = numel( list_nb );

if ~exist( 'list_m', 'var' ) || isempty( list_m ),  list_m = 10 : 5 : 20;  end
if ~exist( 'nrtest', 'var' ) || isempty( nrtest ),  nrtest = 5;  end
nlistm = numel( list_m );

tic

counter  = 0;
tol     = 100*eps^( 2/3 );
tole    = sqrt( eps );  tolm = 10^5*eps;
max_err = 0.0;  failures = 0;

cmpq1 = 1;  orthm = 2;  % Only SVD is used.

if ~exist( 'list_m', 'var' ) || isempty( list_m ),  list_m = 10 : 5 : 20;  end
if ~exist( 'nrtest', 'var' ) || isempty( nrtest ),  nrtest = 5;  end

unsym    = zeros( nrtest, 3 );
unsym_nb = zeros( nrtest, nrnb );
err      = zeros( nrtest, nrnb );
ert      = zeros( nrtest, nrnb );
timing   = zeros( nrtest, nrnb+1 );
timingS  = zeros( nrtest, nrnb+1 );
mtiming  = zeros( nlistm, nrnb+1 );
mtimingS = zeros( nlistm, nrnb+1 );
meanerr  = zeros( nlistm, nrnb );
meanerrS = zeros( nlistm, nrnb );
meanuns  = zeros( nlistm, 3 );
%
% Using rand.
%
for m = list_m,
    n = 2*m;
    counter = counter + 1;
    disp( [ 'n = ', num2str( n ), ',  m = ', num2str( m ) ] )
    eig2 = zeros( n, 1 );
    Ev2  = zeros( n, 1 );
    for j = 1 : nrtest,
        %
        % Generate the input matrices.
        %
        A = 10*rand( m ) - 5;
        D = 10*rand( m ) - 5;
        E = 10*rand( m ) - 5;
        D = D - D';  D = D/2;  E = E - E';
        S = [ A D; E A' ];
        %
        DE = [ tril( E ) zeros( m, 1 ) ] +  [ zeros( m, 1 ) triu( D ) ];
        %
        B = 10*rand( m ) - 5;
        F = 10*rand( m ) - 5;
        G = 10*rand( m ) - 5;
        F = F + F';  F = F/2;
        G = G + G';
        H = [ B F; G -B' ];
        %
        FG = [ tril( G ) zeros( m, 1 ) ] +  [ zeros( m, 1 ) triu( F ) ];
        clear D E F G
        %
        % MATLAB computation.
        %
        time = tic;
        eig1 = eig( H, S );
        time = toc( time );
        timing(j,1) = time;
        %
        % For computing the stable deflating subspace.
        %
        time = tic;
        [ AA, BB, Qq, Z ] = qz( H, S, 'real' );
        Ev = ordeig( AA, BB );
        ii = 1;
        while ii <= n,
           if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
              Ev(ii) = imag( Ev(ii) )*1i;  Ev(ii+1) = -Ev(ii);  ii = ii + 2;
           else
              ii = ii + 1;
           end
        end
        SELECT = real( Ev < 0 );  ns = numel( find( SELECT == 1 ) );
        [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Z, SELECT );
        Ev1  = eig( AAS, BBS );
        time = toc( time );
        timingS(j,1) = time;
        clear AA AAS BB BBS H Qq Qs S SELECT Z
        unsym(j,1) = norm( eig1 - cmpoles( eig1, -eig1 ) );
        unsym(j,2) = norm( Ev1  - cmpoles( Ev1,  -Ev1  ) );
        %
        % FORTRAN computation.
        %
        k = 0;
        for nb = list_nb,
            k = k + 1;
            time = tic;
            [ Alphar, Alphai, Beta ] = skewHamildefl_nb( A, DE, B, FG, nb );
            eig2(1:m) = ( Alphar + 1i*Alphai )./Beta;
            time = toc( time );
            timing(j,1+k) = time;
            eig2(m+1:n) = -eig2(1:m);
            %
            time = tic;
            [ Alphar, Alphai, Beta, Q, neig ] = skewHamildefl_nb( A, DE, B, FG, nb, cmpq1 );
            Ev2(1:m) = ( Alphar + 1i*Alphai )./Beta;
            time = toc( time );
            timingS(j,1+k) = time;
            Ev2(m+1:n) = -Ev2(1:m);
            %
            % Test symmetry property.
            %
            unsym_nb(j,k) = norm( eig2 - cmpoles( eig2, -eig2 ) );
            %
            % Test eigenvalues.
            %
            eivc  = cmpoles( eig1,  eig2 );
            err(j,k) = norm( eig1 - eivc )/max( 1, norm( eig1 ) );
            %
            Evc   = cmpoles( Ev1,  Ev2 );
            ert(j,k) = norm( Ev1 - Evc )/max( 1, norm( Ev1 ) );
            %
            max_err = max( max_err, max( err(j,k), ert(j,k) ) ); 
            %
            Vx = [ Q Zs(:,1:ns) ];  tolr = tolm*n*max( norm( Vx ), 1 );
            ers = max( [ abs( neig - ns ) abs( ns - rank( Vx, tolr ) ) ] );
            if ers > 0,  
               disp( 'Subspace comparison failed 5' )
               failures = failures + 1;
            end
        end
        unsym(:,3)  = mean( unsym_nb, 2 )';  % Mean over nb values (normally, 0).
    end
    if nrtest > 1,
        mtiming (counter,:) = mean( timing  );
        mtimingS(counter,:) = mean( timingS );
        meanuns (counter,:) = mean( unsym   );
        meanerr (counter,:) = mean( err     );
        meanerrS(counter,:) = mean( ert     );
    else
        mtiming (counter,:) = timing;
        mtimingS(counter,:) = timingS;
        meanuns (counter,:) = unsym;
        meanerr (counter,:) = err;
        meanerrS(counter,:) = ert;
    end
    %
    disp( ' ' )
    disp( 'Mean relative errors for eigenvalues and deflating subspace computation with DGHUDP' )
    disp( ' ' )
    disp( '------------------------------------------' )
    disp( '   nb    eigenvalues   deflating subspace ' )
    disp( '------------------------------------------' )
    for jj = 1 : nrnb,
        disp( [ sprintf( '%5d', list_nb(jj) ), ...
                sprintf( '%15.5g', meanerr (counter,jj) ), ...
                sprintf( '%15.5g', meanerrS(counter,jj) ) ] )
    end
    disp( '------------------------------------------' )
    disp( ' ' )
    %
    disp( 'Timings (sec) for eigenvalues and deflating subspace computation' )
    disp( ' ' )
    disp( '---------------------------------' )
    disp( '            eig      eig+ordeig ' )
    disp( '---------------------------------' )
    disp( [ '     ', sprintf( '%13.5g', mtiming (counter,1) ), ...
                     sprintf( '%13.5g', mtimingS(counter,1) ) ] )
    disp( '---------------------------------' )
    disp( '   nb   DGHUDP-eig   DGHUDP-defl ' )
    disp( '---------------------------------' )
    for jj = 1 : nrnb,
        disp( [ sprintf( '%5d', list_nb(jj) ), ...
                sprintf( '%13.5g', mtiming (counter,jj+1) ), ...
                sprintf( '%13.5g', mtimingS(counter,jj+1) ) ] )
    end
    disp( '---------------------------------' )
    disp( ' ' )
end

toc

disp( [ 'Sum of timings = ', sprintf( '%13.5g', nrtest*sum( sum( mtiming + mtimingS ) ) ) ] )
disp( ' ' )

countT = counter*nrtest;
if max_err < tole,
   disp( [ 'DGHUDP :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( countT  ) ] )
   disp( ' ' )
else
   disp( [ 'DGHUDP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( countT  ) ] )
end

disp( 'Mean asymmetry errors for eigenvalue computation' )
disp( ' ' )
disp( '---------------------------------------------------' )
disp( 'Test #   m                 asymmetry error         ' )
disp( '                     eig       ordeig       DGHUDP ' )
disp( '---------------------------------------------------' )
for ii = 1 : counter,
    disp( [ sprintf( '%4d', ii ), sprintf( '%7d', list_m(ii) ), ...
            sprintf( '%13.5g', meanuns(ii,:) ) ] )
end
disp( '---------------------------------------------------' )
disp( ' ' )

currdate = date;
filnam = [ 'perf_skewHamildefl_nb_' currdate '_large' ];
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
save( filnam, 'counter', 'list_m', 'list_nb', 'meanerr', 'meanerrS', 'meanuns', 'mtiming', 'mtimingS', 'nrtest' )
%
disp( [ 'Main results saved in the file ', filnam ] )
