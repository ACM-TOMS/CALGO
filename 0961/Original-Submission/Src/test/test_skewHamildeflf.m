function test_skewHamildeflf
% Function for testing the routine DGHFDF. Random matrices of all
% orders m = 1:50 are used and the relative errors are checked.

%
%   Contributor:
%   V. Sima, Feb. 2011.
%
%   Revisions:
%   M. Voigt, July 2013.
%   V. Sima, July 2014.
%

disp( ' ' )
disp( 'Random tests of the routine DGHFDF' )
disp( ' ' )

count = 0;
max_err = 0.0;
tole = sqrt( eps );  tolm = 10^5*eps;
tol  = 10*tole;

cmpq0 = 0;  cmpu0 = 0;  cmpq1 = 1;  cmpu1 = 1;
orth0 = 0;  orth1 = 1;  orth2 = 2;

% Check default values.

m = 4;  n = 2*m;
%
Z = rand( n );  B = rand( m );  FG = rand( m, m+1 );
%
F = triu( FG(:,2:end) );     G = tril( FG(:,1:m) );
%
J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
S = J*Z'*J'*Z;
H = [ B F+triu( F, 1 )'; G+tril( G, -1 )' -B' ];
%
[ Alphar,   Alphai,   Beta,  Q,  U,  neig  ] = skewHamildeflf( Z, B, FG, cmpq1, cmpu1 );
[ Alphar1,  Alphai1,  Beta1, Q1, U1        ] = skewHamildeflf( Z, B, FG, cmpq1, cmpu1 );
[ Alphar2,  Alphai2,  Beta2, Q2,     neig1 ] = skewHamildeflf( Z, B, FG, cmpq1, cmpu0 );
[ Alphar3,  Alphai3,  Beta3, Q3            ] = skewHamildeflf( Z, B, FG, cmpq1, cmpu0 );
[ Alphar4,  Alphai4,  Beta4, Q4,     neig2 ] = skewHamildeflf( Z, B, FG, cmpq1 );
[ Alphar5,  Alphai5,  Beta5, Q5            ] = skewHamildeflf( Z, B, FG, cmpq1 );
[ Alphar6,  Alphai6,  Beta6,     U2, neig3 ] = skewHamildeflf( Z, B, FG, cmpq0, cmpu1 );
[ Alphar7,  Alphai7,  Beta7,     U3        ] = skewHamildeflf( Z, B, FG, cmpq0, cmpu1 );
[ Alphar8,  Alphai8,  Beta8                ] = skewHamildeflf( Z, B, FG, cmpq0, cmpu0 );
[ Alphar9,  Alphai9,  Beta9                ] = skewHamildeflf( Z, B, FG, cmpq0 );
[ Alphara,  Alphaia,  Betaa                ] = skewHamildeflf( Z, B, FG );
[ Alpharb,  Alphaib,  Betab, Q6, U4, neig4 ] = skewHamildeflf( Z, B, FG, cmpq1, ...
                                                                         cmpu1, orth0 );
[ Alpharc,  Alphaic,  Betac, Q7, U5        ] = skewHamildeflf( Z, B, FG, cmpq1, ...
                                                                         cmpu1, orth0 );
[ Alphard,  Alphaid,  Betad, Q8, U6, neig5 ] = skewHamildeflf( Z, B, FG, cmpq1, ...
                                                                         cmpu1, orth1 );
[ Alphare,  Alphaie,  Betae, Q9, U7        ] = skewHamildeflf( Z, B, FG, cmpq1, ...
                                                                         cmpu1, orth1 );
[ Alpharf,  Alphaif,  Betaf, Qa, U8, neig6 ] = skewHamildeflf( Z, B, FG, cmpq1, ...
                                                                         cmpu1, orth2 );
[ Alpharg,  Alphaig,  Betag, Qb, U9        ] = skewHamildeflf( Z, B, FG, cmpq1, ...
                                                                         cmpu1, orth2 );
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
             norm( Beta   - Beta8    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphar9  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphai9  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Beta9    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphara  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaia  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betaa    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alpharb  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaib  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betab    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alpharc  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaic  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betac    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphard  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaid  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betad    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alphare  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaie  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betae    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alpharf  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaif  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betaf    )/max( 1, norm( Beta   ) ), ...
             norm( Alphar - Alpharg  )/max( 1, norm( Alphar ) ), ...
             norm( Alphai - Alphaig  )/max( 1, norm( Alphai ) ), ...
             norm( Beta   - Betag    )/max( 1, norm( Beta   ) ) ] );
err = max( [ neig - neig1, neig - neig2, neig - neig3, neig - neig4, ...
             neig - neig5, neig - neig6, err ] );
%
erq = max( [ norm( Q  - Q1 ), norm( Q - Q2 ), norm( Q  - Q3 ), norm( Q  - Q4 ), ...
             norm( Q  - Q5 ), norm( Q - Q6 ), norm( Q  - Q7 ), norm( Q8 - Q9 ), ... 
             norm( Qa - Qb ), norm( U - U1 ), norm( U  - U2 ), norm( U  - U3 ), ... 
             norm( U  - U4 ), norm( U - U5 ), norm( U6 - U7 ), norm( U8 - U9 ) ] )/n;
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
for k = 1 : 5,
   for m = 0 : 50,
      n = 2*m;
      Z = rand( n );  B = rand( m );  FG = rand( m, m+1 );
      %
      F = triu( FG(:,2:end) );     G = tril( FG(:,1:m) );
      %
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      S = J*Z'*J'*Z;
      H = [ B F+triu( F, 1 )'; G+tril( G, -1 )' -B' ];
      %
      evm = eig( H, S );
      if isreal( evm ),
         evp = evm( evm > 0 );
      else
         evp = union( evm( imag( evm ) == 0 & real( evm >= 0 ) ), ...
                      evm( imag( evm ) > 0 ) );
      end
      %
      % For checking the deflating subspace.
      %
      [ AA, BB, Qq, Zz ] = qz( H, S, 'real' );
      Ev = ordeig( AA, BB );
      ii = 1;
      while ii <= n,
         if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
            Ev(ii) = imag( Ev(ii) )*1i;  Ev(ii+1) = -Ev(ii);  ii = ii + 2;
         else
            ii = ii + 1;
         end
      end
      SELECT = real( Ev < 0 );  ns = sum( SELECT );
      [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Zz, SELECT );
      %
      count = count + 1;
      for orthm = 0 : 2,
         [ Alphar, Alphai, Beta, Q, U, neig ] = skewHamildeflf( Z, B, FG, ...
                                                    cmpq1, cmpu1, orthm );
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
   disp( [ 'DGHFDF :    passed  -- maximum relative error max_err = ', ...
                        num2str( max_err ) ] )
   disp( [ '            Number of problems solved using rand      = ', ...
                        num2str( count   ) ] )
   disp( ' ' )
else
   disp( [ 'DGHFDF :    failed  -- maximum relative error max_err = ', ...
                        num2str( max_err ) ] )
   disp( [ '            Number of problems solved using rand      = ', ...
                        num2str( count   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', ...
                        num2str( failures ) ] )
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
for k = 1 : 5,
   for m = 0 : 50,
      n = 2*m;
      Z = randn( n );  B = randn( m );  FG = randn( m, m+1 );
      %
      F = triu( FG(:,2:end) );     G = tril( FG(:,1:m) );
      %
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      S = J*Z'*J'*Z;
      H = [ B F+triu( F, 1 )'; G+tril( G, -1 )' -B' ];
      %
      evm = eig( H, S );
      if isreal( evm ),
         evp = evm( evm > 0 );
      else
         evp = union( evm( imag( evm ) == 0 & real( evm >= 0 ) ), ...
                      evm( imag( evm ) > 0 ) );
      end
      %
      % For checking the deflating subspace.
      %
      [ AA, BB, Qq, Zz ] = qz( H, S, 'real' );
      Ev = ordeig( AA, BB );
      ii = 1;
      while ii <= n,
         if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
            Ev(ii) = imag( Ev(ii) )*1i;  Ev(ii+1) = -Ev(ii);  ii = ii + 2;
         else
            ii = ii + 1;
         end
      end
      SELECT = real( Ev < 0 );  ns = sum( SELECT );
      [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Zz, SELECT );
      %
      count = count + 1;
      for orthm = 0 : 2,
         [ Alphar, Alphai, Beta, Q, U, neig ] = skewHamildeflf( Z, B, FG, ...
                                                    cmpq1, cmpu1, orthm );
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
   disp( [ 'DGHFDF :    passed  -- maximum relative error max_err = ', ...
                        num2str( max_err ) ] )
   disp( [ '            Number of problems solved using randn     = ', ...
                        num2str( count   ) ] )
   disp( ' ' )
else
   disp( [ 'DGHFDF :    failed  -- maximum relative error max_err = ', ...
                        num2str( max_err ) ] )
   disp( [ '            Number of problems solved using randn     = ', ...
                        num2str( count   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', ...
                        num2str( failures ) ] )
   disp( ' ' )
end        
