function test_skewHamil2eig_nb
% Function for testing the routine DGHUSP. Random matrices of all
% orders m = 1:20 are used and the relative errors are checked.

%
%   Contributor:
%   V. Sima, Dec. 2010.
%
%   Revisions:
%   V. Sima, Jul. 2013.
%   M. Voigt, Jul. 2013.
%

disp( ' ' )
disp( 'Random tests of the routine DGHUSP' )
disp( ' ' )

if ~exist( 'nb', 'var' ) || isempty( nb ),  nb = 4;  end

count   = 0;
max_err = 0.0;
tol  = 100*eps^( 2/3 );
tole = sqrt( eps );

job0  = 0;  job1  = 1;
cmpq0 = 0;  cmpq1 = 1;  cmpq2 = 2;

%
% Check functionality. A large m is also used.
%
% Using rand.
%
for k = 1 : 2,
   for m = [ 0 : 20, 210 ]
      nbmax = min( max( m, 1 ), 20 );
      n = 2*m;
      A = rand( m );  DE = rand( m, m+1 );
      B = rand( m );  FG = rand( m, m+1 );
      %
      D = triu( DE(:,2:end), 1 );  E = tril( DE(:,1:m), -1 );
      F = triu( FG(:,2:end), 1 );  G = tril( FG(:,1:m), -1 );
      %
      S = [ A D-D'; E-E' A' ];
      T = [ B F-F'; G-G' B' ];
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      %
      evm = eig( T, S );
      %
      for job = 0 : 1,
         count = count + 1;
         [ Alphar, Alphai, Beta, Ao, Do, Bo, Fo, Q ] = ...
                            skewHamil2eig_nb( A, DE, B, FG, nb, job, cmpq1 );
         ev  = ( Alphar + Alphai*1i )./Beta;
         evd = [ ev; ev ];
         err = norm( ev - cmpoles( ev, evd ) )/max( 1, norm( ev ) );
         if err > tol,  disp( 'Failed 2' ),  return,  end
         max_err = max( max_err, err );
         if job == 1,
            Dot = triu( Do, 1 );
            Fot = triu( Fo, 1 );
            %
            So = [ Ao Dot-Dot'; zeros( m ) Ao' ];
            To = [ Bo Fot-Fot'; zeros( m ) Bo' ]; 
            errm = max( [ norm( So - J*Q'*J'*S*Q, 1 )/max( norm( So, 1 ), 1 ), ...
                          norm( To - J*Q'*J'*T*Q, 1 )/max( norm( To, 1 ), 1 ) ] );
            if errm > tol,  disp( 'Failed 3' ),  return,  end
            max_err = max( max_err, errm );
         end
      end
   end
end
%
% Using randn.
%
for k = 1 : 2,
   for m = [ 0 : 20, 210 ]
      nbmax = min( max( m, 1 ), 20 );
      n = 2*m;
      A = randn( m );  DE = randn( m, m+1 );
      B = randn( m );  FG = randn( m, m+1 );
      %
      D = triu( DE(:,2:end), 1 );  E = tril( DE(:,1:m), -1 );
      F = triu( FG(:,2:end), 1 );  G = tril( FG(:,1:m), -1 );
      %
      S = [ A D-D'; E-E' A' ];
      T = [ B F-F'; G-G' B' ];
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      %
      evm = eig( T, S );
      %
      for job = 0 : 1,
         count = count + 1;
         [ Alphar, Alphai, Beta, Ao, Do, Bo, Fo, Q ] = ...
                            skewHamil2eig_nb( A, DE, B, FG, nb, job, cmpq1 );
         ev  = ( Alphar + Alphai*1i )./Beta;
         evd = [ ev; ev ];
         err = norm( ev - cmpoles( ev, evd ) )/max( 1, norm( ev ) );
         if err > tol,  disp( 'Failed 4' ),  return,  end
         max_err = max( max_err, err );
         if job == 1,
            Dot = triu( Do, 1 );
            Fot = triu( Fo, 1 );
            %
            So = [ Ao Dot-Dot'; zeros( m ) Ao' ];
            To = [ Bo Fot-Fot'; zeros( m ) Bo' ];
            errm = max( [ norm( So - J*Q'*J'*S*Q, 1 )/max( norm( So, 1 ), 1 ), ...
                          norm( To - J*Q'*J'*T*Q, 1 )/max( norm( To, 1 ), 1 ) ] );
            if errm > tol,  disp( 'Failed 5' ),  return,  end
            max_err = max( max_err, errm );
         end
      end
   end
end

if max_err < tole,
   disp( [ 'DGHUSP :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( count   ) ] )
   disp( ' ' )
else
   disp( [ 'DGHUSP :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( count   ) ] )
end
