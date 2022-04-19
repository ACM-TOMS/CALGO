function test_skewHamil2eig
% Function for testing the routine DGHUST. Random matrices of all
% orders m = 1:20 are used and the relative errors are checked.

%
%   Contributor:
%   V. Sima, Dec. 2010.
%
%   Revisions:
%   M. Voigt, Jul. 2013.
%

disp( ' ' )
disp( 'Random tests of the routine DGHUST' )
disp( ' ' )

count   = 0;
max_err = 0.0;
tol  = 100*eps^( 2/3 );
tole = sqrt( eps );

job0  = 0;  job1  = 1;
cmpq0 = 0;  cmpq1 = 1;  cmpq2 = 2;

% Check default values.

m = 4;  n = 2*m;
%
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
[ Q, r ] = qr( S );
%
[ Alphar,   Alphai,   Beta,   Ao,   Do,   Bo,  Fo,  Qr ] = ...
                            skewHamil2eig( A, DE, B, FG, job1, cmpq2, Q );
[ Alphar1,  Alphai1,  Beta1,  Ao1,  Do1,  Bo1, Fo1, Q1 ] = ...
                            skewHamil2eig( A, DE, B, FG, job1, cmpq1 );
[ Alphar2,  Alphai2,  Beta2,  Ao2,  Do2,  Bo2, Fo2 ] = ...
                            skewHamil2eig( A, DE, B, FG, job1, cmpq0 );
[ Alphar3,  Alphai3,  Beta3,  Ao3,  Do3,  Bo3, Fo3 ] = ...
                            skewHamil2eig( A, DE, B, FG, job1 );
[ Alphar4,  Alphai4,  Beta4,  Ao4,  Do4,  Bo4 ] = ...
                            skewHamil2eig( A, DE, B, FG, job1 );
[ Alphar5,  Alphai5,  Beta5,  Ao5,  Do5 ] = ...
                            skewHamil2eig( A, DE, B, FG, job1 );
[ Alphar6,  Alphai6,  Beta6,  Ao6 ] = ...
                            skewHamil2eig( A, DE, B, FG, job1 );
[ Alphar7,  Alphai7,  Beta7 ] = ...
                            skewHamil2eig( A, DE, B, FG, job0 );
[ Alphar8 , Alphai8,  Beta8 ] = ...
                            skewHamil2eig( A, DE, B, FG );
%
Do  = triu( Do,  1 );  Do1 = triu( Do1, 1 );  Do2 = triu( Do2, 1 );
Do3 = triu( Do3, 1 );  Do4 = triu( Do4, 1 );  Do5 = triu( Do5, 1 );
Fo  = triu( Fo,  1 );  Fo1 = triu( Fo1, 1 );  Fo2 = triu( Fo2, 1 );
Fo3 = triu( Fo3, 1 );
%
err = max( [ norm( Alphar - Alphar1  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai1  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta1    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar2  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai2  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta2    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar3  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai3  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta3    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar4  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai4  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta4    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar5  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai5  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta5    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar6  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai6  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta6    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar7  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai7  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta7    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar8  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai8  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta8    )/max( 1, norm( Beta   ) ) ] );
%
err = max( [ norm( Ao  - Ao1,  1 )/max( 1, norm( Ao,  1 ) ), ...
             norm( Do  - Do1,  1 )/max( 1, norm( Do,  1 ) ), ...
             norm( Bo  - Bo1,  1 )/max( 1, norm( Bo,  1 ) ), ...
             norm( Fo  - Fo1,  1 )/max( 1, norm( Fo,  1 ) ), ...
             norm( Ao  - Ao2,  1 )/max( 1, norm( Ao,  1 ) ), ...
             norm( Do  - Do2,  1 )/max( 1, norm( Do,  1 ) ), ...
             norm( Bo  - Bo2,  1 )/max( 1, norm( Bo,  1 ) ), ...
             norm( Fo  - Fo2,  1 )/max( 1, norm( Fo,  1 ) ), ...
             norm( Ao  - Ao3,  1 )/max( 1, norm( Ao,  1 ) ), ...
             norm( Do  - Do3,  1 )/max( 1, norm( Do,  1 ) ), ...
             norm( Bo  - Bo3,  1 )/max( 1, norm( Bo,  1 ) ), ...
             norm( Fo  - Fo3,  1 )/max( 1, norm( Fo,  1 ) ), ...
             norm( Ao  - Ao4,  1 )/max( 1, norm( Ao,  1 ) ), ...
             norm( Do  - Do4,  1 )/max( 1, norm( Do,  1 ) ), ...
             norm( Bo  - Bo4,  1 )/max( 1, norm( Bo,  1 ) ), ...
             norm( Ao  - Ao5,  1 )/max( 1, norm( Ao,  1 ) ), ...
             norm( Do  - Do5,  1 )/max( 1, norm( Do,  1 ) ), ...
             norm( Ao  - Ao6,  1 )/max( 1, norm( Ao,  1 ) ), ...
             err ] );
erq = norm( Q*Q1 - Qr )/n;
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
% Check functionality.
%
% Using rand.
%
for k = 1 : 100,
   for m = 0 : 20,
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
                            skewHamil2eig( A, DE, B, FG, job, cmpq1 );
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
for k = 1 : 100,
   for m = 0 : 20,
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
                            skewHamil2eig( A, DE, B, FG, job, cmpq1 );
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
   disp( [ 'DGHUST :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( count   ) ] )
   disp( ' ' )
else
   disp( [ 'DGHUST :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved                 = ', num2str( count   ) ] )
end
