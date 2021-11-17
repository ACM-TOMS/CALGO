function perf_skewHamil2eig_nb
% Function for testing the performance of the routine DGHUSP (block code).
% Random matrices of orders 2*list_m are used and the relative errors are checked.
% For each m in list_m, nrtest are performed, and the mean values are
% returned.  Default values: list_m = 10 : 5 : 20;  nrtest = 5;

%
%   Contributor:
%   V. Sima, Jul. 2013.
%
%   Revisions:
%   M. Voigt, Jul. 2013.
%

disp( ' ' )
disp( 'Random performance tests of the routine DGHUSP' )
disp( ' ' )

if ~exist( 'list_nb', 'var' ) || isempty( list_nb ),  list_nb = 32;  end
nrnb = numel( list_nb );

if ~exist( 'list_m', 'var' ) || isempty( list_m ),  list_m = 10 : 5 : 20;  end
if ~exist( 'nrtest', 'var' ) || isempty( nrtest ),  nrtest = 5;  end

tic

tol   = 100*eps^( 2/3 );
tole  = sqrt( eps );  tolm = 10^5*eps;
max_err = 0.0;  failures = 0;

job = 1;  cmpq1 = 1;

unsym    = zeros( nrtest, 2 );
err      = zeros( nrtest, nrnb );
ert      = zeros( nrtest, nrnb );
timing   = zeros( nrtest, nrnb+1 );
timingS  = zeros( nrtest, nrnb+1 );  nlistm = numel( list_m );
mtiming  = zeros( nlistm, nrnb+1 );
mtimingS = zeros( nlistm, nrnb+1 );
meanerr  = zeros( nlistm, nrnb );
meanerrS = zeros( nlistm, nrnb );
meanuns  = zeros( nlistm, 2 );
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
        A = rand( m );  DE = rand( m, m+1 );
        B = rand( m );  FG = rand( m, m+1 );
        %
        D = triu( DE(:,2:end), 1 );  E = tril( DE(:,1:m), -1 );
        F = triu( FG(:,2:end), 1 );  G = tril( FG(:,1:m), -1 );
        %
        S = [ A D-D'; E-E' A' ];
        H = [ B F-F'; G-G' B' ];
        clear D E F G
        %
        % MATLAB computation.
        %
        time = tic;
        eig1 = eig( H, S );
        time = toc( time );
        timing(j,1) = time;
        %
        time = tic;
        [ AA, BB, Qq, Z ] = qz( H, S, 'real' );
        Ev1 = ordeig( AA, BB );
        time = toc( time );
        timingS(j,1) = time;
        clear AA BB H Qq S Z
        %
        % FORTRAN computation.
        %
        k = 0;
        for nb = list_nb,
            k = k + 1;
            time = tic;
            [ Alphar, Alphai, Beta ] = skewHamil2eig_nb( A, DE, B, FG, nb );
            eig2 = ( Alphar + 1i*Alphai )./Beta;
            time = toc( time );
            timing(j,1+k) = time;
            eig2 = [ eig2; eig2 ];
            %
            time = tic;
            [ Alphar, Alphai, Beta, Ao, Do, Bo, Fo, Q ] = ...
                      skewHamil2eig_nb( A, DE, B, FG, nb, job, cmpq1 );
            Ev2  = ( Alphar + 1i*Alphai )./Beta;
            time = toc( time );
            timingS(j,1+k) = time;
            Ev2 = [ Ev2; Ev2 ];
            %
            % Test eigenvalues.
            %
            eivc  = cmpoles( eig2,  eig1 );
            err(j,k) = norm( eig2 - eivc )/max( 1, norm( eig2 ) );
            %
            Evc   = cmpoles( Ev2,  Ev1 );
            ert(j,k) = norm( Ev2 - Evc )/max( 1, norm( Ev2 ) );
            %
            max_err = max( max_err, max( err(j,k), ert(j,k) ) ); 
            %
            % Test the multiplicity for eig (it should be doubled).
            %
            eig1s = cmpoles( eivc, eig1 );
            unsym(j,1) = norm( eig1s(1:m) - eig1s(m+1:n) )/norm( eig1s(1:m) );
            eig2s = cmpoles( Evc, Ev1 );
            unsym(j,2) = norm( eig2s(1:m) - eig2s(m+1:n) )/norm( eig2s(1:m) );
        end
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
    disp( 'Mean relative errors for eigenvalue and transformations computation with DGHUSP' )
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
    disp( '   nb             DGHUSP                  ' )
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
   disp( [ 'DGHUSP :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( countT  ) ] )
   disp( ' ' )
else
   disp( [ 'DGHUSP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( countT  ) ] )
end
    
disp( 'Mean ''asymmetry'' errors for eigenvalue computation' )
disp( ' ' )
disp( '--------------------------------------' )
disp( 'Test #   m           asymmetry error  ' )
disp( '                     eig       ordeig ' )
disp( '--------------------------------------' )
for ii = 1 : counter,
    disp( [ sprintf( '%4d', ii ), sprintf( '%7d', list_m(ii) ), ...
        sprintf( '%13.5g', meanuns(ii,:) ) ] )
end
disp( '--------------------------------------' )
disp( ' ' )

currdate = date;
filnam = [ 'perf_skewHamil2eig_' currdate '_large' ];
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
save( filnam, 'counter', 'list_m', 'list_nb', 'meanerr', 'meanerrS', 'mtiming', 'mtimingS', 'nrtest' )
%
disp( [ 'Main results saved in the file ', filnam ] )
