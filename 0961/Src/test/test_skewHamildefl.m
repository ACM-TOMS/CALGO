function test_skewHamildefl
% Function for testing the routine DGHUDF. Random matrices of all
% orders m = 1:100 are used and the relative errors are checked.

%
%   Contributor:
%   V. Sima, Oct. 2009.
%
%   Revisions:
%   V. Sima, Nov. 2010, Dec. 2010, June 2014.
%   M. Voigt, July 2013.
%

disp(' ')
disp('Random tests of the routine DGHUDF')
disp(' ')

count = 0;
max_err = 0.0;
tole = sqrt( eps );  tolm = 10^5*eps;
tol  = 10*tole;

cmpq0 = 0;  cmpq1 = 1;
orth0 = 0;  orth1 = 1;  orth2 = 2;

% Check default values.

m = 4;  n = 2*m;
%
A = rand( m );  DE = rand( m, m+1 );
B = rand( m );  FG = rand( m, m+1 );
%
D = triu( DE(:,2:end), 1 );  E = tril( DE(:,1:m), -1 );
F = triu( FG(:,2:end) );     G = tril( FG(:,1:m) );
%
S = [ A D-D'; E-E' A' ];
H = [ B F+triu( F, 1 )'; G+tril( G, -1 )' -B' ];
J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
%
[ Alphar,   Alphai,   Beta,  Q,  neig  ] = skewHamildefl( A, DE, B, FG, cmpq1 );
[ Alphar1,  Alphai1,  Beta1, Q1        ] = skewHamildefl( A, DE, B, FG, cmpq1 );
[ Alphar2,  Alphai2,  Beta2            ] = skewHamildefl( A, DE, B, FG, cmpq0 );
[ Alphar3,  Alphai3,  Beta3            ] = skewHamildefl( A, DE, B, FG );
[ Alphar4,  Alphai4,  Beta4, Q2, neig1 ] = skewHamildefl( A, DE, B, FG, cmpq1, orth0 );
[ Alphar5,  Alphai5,  Beta5, Q3, neig2 ] = skewHamildefl( A, DE, B, FG, cmpq1, orth1 );
[ Alphar6,  Alphai6,  Beta6, Q4, neig3 ] = skewHamildefl( A, DE, B, FG, cmpq1, orth2 );
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
             norm( Beta   - Beta6    )/max( 1, norm( Beta   ) ) ] );
err = max( [ neig - neig1, neig - neig2, neig - neig3, err ] );
%
erq = norm( Q - Q1 )/n;
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
failures = 0;
%
for k = 1 : 100,
   for m = 0 : 20,
      n = 2*m;
      A = rand( m );  DE = rand( m, m+1 );
      B = rand( m );  FG = rand( m, m+1 );
      %
      D = triu( DE(:,2:end), 1 );  E = tril( DE(:,1:m), -1 );
      F = triu( FG(:,2:end) );     G = tril( FG(:,1:m) );
      %
      S = [ A D-D'; E-E' A' ];
      H = [ B F+triu( F, 1 )'; G+tril( G, -1 )' -B' ];
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      %
      evm = eig( H, S );
      if isreal( evm ),
         evp = evm( evm > 0 );
      else
         evp = union( evm( imag( evm ) == 0 & real( evm > 0 ) ), ...
                      evm( imag( evm ) > 0 ) );
      end
      %
      % For checking the deflating subspace.
      %
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
      %
      count = count + 1;
      for orthm = 0 : 2,
         [ Alphar, Alphai, Beta, Q, neig ] = skewHamildefl( A, DE, B, FG, cmpq1, orthm );
         ev  = ( Alphar + Alphai*1i )./Beta;
         err = norm( ev - cmpoles( ev, evp ) )/max( 1, norm( evp ) );
         if err > tol,  
            disp( [ 'Failed 2, orthm = ', num2str( orthm ) ] )
            failures = failures + 1;
         else
            max_err = max( max_err, err );
         end
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
if max_err < tol,
   disp( [ 'DGHUDF :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using rand      = ', num2str( count   ) ] )
   disp( ' ' )
else
   disp( [ 'DGHUDF :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using rand      = ', num2str( count   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
   disp( ' ' )
end        
%
count = 0;
max_err = 0.0;
%
% Using randn.
%
failures = 0;
%
for k = 1 : 100,
   for m = 0 : 20,
      n = 2*m;
      A = randn( m );  DE = randn( m, m+1 );
      B = randn( m );  FG = randn( m, m+1 );
      %
      D = triu( DE(:,2:end), 1 );  E = tril( DE(:,1:m), -1 );
      F = triu( FG(:,2:end) );     G = tril( FG(:,1:m) );
      %
      S = [ A D-D'; E-E' A' ];
      H = [ B F+triu( F, 1 )'; G+tril( G, -1 )' -B' ];
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      %
      evm = eig( H, S );
      if isreal( evm ),
         evp = evm( evm > 0 );
      else
         evp = union( evm( imag( evm ) == 0 & real( evm > 0 ) ), ...
                      evm( imag( evm ) > 0 ) );
      end
      %
      % For checking the deflating subspace.
      %
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
      %
      count = count + 1;
      for orthm = 0 : 2,
         [ Alphar, Alphai, Beta, Q, neig ] = skewHamildefl( A, DE, B, FG, cmpq1, orthm );
         ev  = ( Alphar + Alphai*1i )./Beta;
         err = norm( ev - cmpoles( ev, evp ) )/max( 1, norm( evp ) );
         if err > tol,  
            disp( [ 'Failed 4, orthm = ', num2str( orthm ) ] )
            failures = failures + 1;
         else
            max_err = max( max_err, err );
         end
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
if max_err < tol,
   disp( [ 'DGHUDF :    passed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using randn     = ', num2str( count   ) ] )
   disp( ' ' )

else
   disp( [ 'DGHUDF :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using randn     = ', num2str( count   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
   disp( ' ' )
end        
