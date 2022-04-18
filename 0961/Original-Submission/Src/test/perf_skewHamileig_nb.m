function perf_skewHamileig_nb
% Function for testing the performance of the routine DGHUTP (block code).
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
disp( 'Random performance tests of the routine DGHUTP' )
disp( ' ' )

if ~exist( 'list_nb', 'var' ) || isempty( list_nb ),  list_nb = 32;  end
nrnb = numel( list_nb );

if ~exist( 'list_m', 'var' ) || isempty( list_m ),  list_m = 10 : 5 : 20;  end
if ~exist( 'nrtest', 'var' ) || isempty( nrtest ),  nrtest = 5;  end

tic

tol   = 100*eps^( 2/3 );
tole  = sqrt( eps );  tolm = 10^5*eps;
max_err = 0.0;  failures = 0;

job = 1;  cmpq1 = 1;  cmpq2 = 1;  % orth0 = 0;  orth1 = 1;  orth2 = 2;

unsym    = zeros( nrtest, 3 );
unsym_nb = zeros( nrtest, nrnb );
err      = zeros( nrtest, nrnb );
ert      = zeros( nrtest, nrnb );
timing   = zeros( nrtest, nrnb+1 );
timingS  = zeros( nrtest, nrnb+1 );  nlistm = numel( list_m );
mtiming  = zeros( nlistm, nrnb+1 );
mtimingS = zeros( nlistm, nrnb+1 );
meanerr  = zeros( nlistm, nrnb );
meanerrS = zeros( nlistm, nrnb );
meanuns  = zeros( nlistm, 3 );
%
counter  = 0;
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
        Ev1 = ordeig( AA, BB );
        time = toc( time );
        timingS(j,1) = time;
        clear AA BB H Qq S Z
        unsym(j,1) = norm( eig1 - cmpoles( eig1, -eig1 ) );
        unsym(j,2) = norm( Ev1  - cmpoles( Ev1,  -Ev1  ) );
        %
        % FORTRAN computation.
        %
        k = 0;
        for nb = list_nb,
            k = k + 1;
            time = tic;
            [ Alphar, Alphai, Beta ] = skewHamileig_nb( A, DE, B, FG, nb );
            eig2 = ( Alphar + 1i*Alphai )./Beta;
            time = toc( time );
            timing(j,1+k) = time;
            eig2 = [ eig2; -eig2 ];
            %
            time = tic;
            [ Alphar, Alphai, Beta, Ao, Do, Bo, Fo, C1o, Vo, C2o, Q1, Q2 ] = ...
                      skewHamileig_nb( A, DE, B, FG, nb, job, cmpq1, cmpq2 );
            Ev2  = ( Alphar + 1i*Alphai )./Beta;
            time = toc( time );
            timingS(j,1+k) = time;
            Ev2 = [ Ev2; -Ev2 ];
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
    disp( 'Mean relative errors for eigenvalue and transformations computation with DGHUTP' )
    disp( ' ' )
    disp( '------------------------------------------' )
    disp( '   nb    eigenvalues   + transformations  ' )
    disp( '------------------------------------------' )
    for jj = 1 : nrnb,
        disp( [ sprintf( '%5d', list_nb(jj) ), ...
                sprintf( '%15.5g', meanerr (counter,jj) ), ...
                sprintf( '%15.5g', meanerrS(counter,jj) ) ] )
    end
    disp( '------------------------------------------' )
    disp( ' ' )
    %
    disp( 'Timings (sec) for eigenvalue computation and transformations computation' )
    disp( ' ' )
    disp( '------------------------------------------' )
    disp( '       eigenvalues     + transformations  ' )
    disp( '------------------------------------------' )
    disp( '            eig        eig+ordeig         ' )
    disp( '------------------------------------------' )
    disp( [ '     ', sprintf( '%13.5g', mtiming (counter,1) ), '  ', ...
                     sprintf( '%13.5g', mtimingS(counter,1) ) ] )
    disp( '------------------------------------------' )
    disp( '   nb             DGHUTP                  ' )
    disp( '------------------------------------------' )
    for jj = 1 : nrnb,
        disp( [ sprintf( '%5d', list_nb(jj) ), ...
                sprintf( '%13.5g', mtiming (counter,jj+1) ), '  ', ...
                sprintf( '%13.5g', mtimingS(counter,jj+1) ) ] )
    end
    disp( '------------------------------------------' )
    disp( ' ' )
end

toc

disp( [ 'Sum of timings = ', sprintf( '%13.5g', nrtest*sum( sum( mtimingS + mtiming ) ) ) ] )

countT = counter*nrtest;
disp( ' ' )
if max_err < tole,
   disp( [ 'DGHUTP :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( countT  ) ] )
   disp( ' ' )
else
   disp( [ 'DGHUTP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( countT  ) ] )
end
    
disp( 'Mean asymmetry errors for eigenvalue computation' )
disp( ' ' )
disp( '---------------------------------------------------' )
disp( 'Test #   m                 asymmetry error         ' )
disp( '                     eig       ordeig       DGHUTP ' )
disp( '---------------------------------------------------' )
for ii = 1 : counter,
    disp( [ sprintf( '%4d', ii ), sprintf( '%7d', list_m(ii) ), ...
            sprintf( '%13.5g', meanuns(ii,:) ) ] )
end
disp( '---------------------------------------------------' )
disp( ' ' )

currdate = date;
filnam = [ 'perf_skewHamileig_' currdate '_large' ];
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

