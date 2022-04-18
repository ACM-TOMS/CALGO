function test_skewHamildeflfZ
% Function for testing the routine ZGHFDF. Random matrices of all
% orders m = 1:50 are used and the relative errors are checked.

%
%   Contributor:
%   V. Sima, Dec. 2010.
%
%   Revisions:
%   M. Voigt, July 2013.
%   V. Sima, June 2014, July 2014.
%

disp(' ')
disp('Random tests of the routine ZGHFDF')
disp(' ')

Details = 0;  % Set it to 1 for printing failure messages.

count = 0;
tole  = sqrt( eps );  tolm = 10^5*eps;
tol   = 10*tole;
max_err = 0.0;  failures = 0;

cmpq0 = 0;  cmpq1 = 1;  
cmpu0 = 0;  cmpu1 = 1;
orth0 = 0;  orth1 = 1;  orth2 = 2;

% Check default values.

m = 4;  n = 2*m;
%
Z = rand( n ) + rand( n )*1i;  
B = rand( m ) + rand( m )*1i;  FG = rand( m, m+1 ) + rand( m, m+1 )*1i;
%
[ Alpha,  Beta,  Q,  U,  neig  ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu1, orth0 );
[ Alpha1, Beta1, Q1, U1, neig1 ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu1, orth1 );
[ Alpha2, Beta2, Q2, U2, neig2 ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu1, orth2 );
[ Alpha3, Beta3, Q3, U3, neig3 ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu1 );
[ Alpha4, Beta4, Q4,     neig4 ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu0 );
[ Alpha5, Beta5, Q5,     neig5 ] = skewHamildeflfZ( Z, B, FG, cmpq1 );
[ Alpha6, Beta6                ] = skewHamildeflfZ( Z, B, FG, cmpq0 );
[ Alpha7, Beta7                ] = skewHamildeflfZ( Z, B, FG );
%
err = max( [ norm( Alpha - Alpha1, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha2, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha3, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha4, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha5, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha6, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Alpha - Alpha7, 1 )/max( 1, norm( Alpha, 1 ) ), ...
             norm( Beta  - Beta1,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta2,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta3,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta4,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta5,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta6,  1 )/max( 1, norm( Beta,  1 ) ), ...
             norm( Beta  - Beta7,  1 )/max( 1, norm( Beta,  1 ) ) ] );
%
err = max( [ neig - neig1, neig - neig2, neig - neig3, neig - neig4, neig - neig5, err ] );
%
erq = max( [ norm( Q - Q3 ) norm( Q4 - Q5 ) norm( U - U3 ) ] )/n;
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
for k = 1 : 5,
   for m = 0 : 50,
      n = 2*m;
      Z = rand( n ) + rand( n )*1i; 
      B = rand( m ) + rand( m )*1i;  FG = rand( m, m+1 ) + rand( m, m+1 )*1i;
      v = diag( FG(:,1:m   ) );  FG(1 : m+1:m*m    ) = ( v + conj( v ) )/2;
      v = diag( FG(:,2:m+1 ) );  FG(m+1:m+1:m*(m+1)) = ( v + conj( v ) )/2;
      F = triu( FG(:,2:end) );   F = F + triu( F,  1 )';
      G = tril( FG(:,1:m) );     G = G + tril( G, -1 )';
      %
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      S = J*Z'*J'*Z;
      H = [ B F; G -B' ];
      %
      if m == 0,  S = [ ];  H = [ ];  end
      evm = eig( H, S );
      %
      count = count + 1;
      %
      % For checking the deflating subspace.
      %
      [ AA, BB, Qq, Zz ] = qz( H, S );
      Ev = ordeig( AA, BB );
      for ii = 1 : n,
         if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
            Ev(ii) = imag( Ev(ii) )*1i;
         end
      end
      SELECT = real( Ev < 0 );  ns = sum( SELECT );
      [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Zz, SELECT );
      %
      for orthm = 0 : 2,
         [ Alpha, Beta, Q, U, neig ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu1, orthm );
         eva = Alpha./Beta;  eva = eva(:);
         err = norm( evm - cmpoles( evm, eva ) )/max( 1, norm( evm ) );
         if err > max( 1, n )*tol,
            if Details == 1,  disp( 'Failed 2' ),  end
            failures = failures + 1;
         else
            max_err = max( max_err, err );
         end
         %
         Vx = [ Q Zs(:,1:ns) ];  tolr = tolm*n*max( norm( Vx ), 1 );
         err = max( [ abs( neig - ns ) abs( ns - rank( Vx, tolr ) ) ] );
         if err > 0,  
            if Details == 1,  disp( [ 'Failed 3, orthm = ', num2str( orthm ) ] ),  end
            failures = failures + 1;
         end
      end
   end
end

if max_err < tol,
   disp( [ 'ZGHFDF :    passed  -- maximum relative error max_err = ', num2str( max_err  ) ] )
   disp( [ '            Number of problems solved using rand      = ', num2str( count    ) ] )
   disp( ' ' )
else
   disp( [ 'ZGHFDF :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using rand      = ', num2str( count   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
   disp( ' ' )
end        
%
count   = 0;
max_err = 0.0;  failures = 0;
%
% Using randn.
%
for k = 1 : 5,
   for m = 0 : 50,
      n = 2*m;
      Z = randn( n ) + randn( n )*1i;
      B = randn( m ) + randn( m )*1i;  FG = randn( m, m+1 ) + randn( m, m+1 )*1i;
      v = diag( FG(:,1:m   ) );  FG(1 : m+1:m*m    ) = ( v + conj( v ) )/2;
      v = diag( FG(:,2:m+1 ) );  FG(m+1:m+1:m*(m+1)) = ( v + conj( v ) )/2;
      F = triu( FG(:,2:end) );   F = F + triu( F,  1 )';
      G = tril( FG(:,1:m) );     G = G + tril( G, -1 )';
      %
      J = [ zeros( m ) eye( m ); -eye( m ) zeros( m ) ];
      S = J*Z'*J'*Z;
      H = [ B F; G -B' ];
      %
      if m == 0,  S = [ ];  H = [ ];  end
      evm = eig( H, S );
      %
      count = count + 1;
      %
      % For checking the deflating subspace.
      %
      [ AA, BB, Qq, Zz ] = qz( H, S );
      Ev = ordeig( AA, BB );
      for ii = 1 : n,
         if abs( real( Ev(ii) ) ) / abs( Ev(ii) ) < tolm*n,
            Ev(ii) = imag( Ev(ii) )*1i;
         end
      end
      SELECT = real( Ev < 0 );  ns = sum( SELECT );
      [ AAS, BBS, Qs, Zs ] = ordqz( AA, BB, Qq, Zz, SELECT );
      %
      for orthm = 0 : 2,
         [ Alpha, Beta, Q, U, neig ] = skewHamildeflfZ( Z, B, FG, cmpq1, cmpu1, orthm );
         eva = Alpha./Beta;  eva = eva(:);
         err = norm( evm - cmpoles( evm, eva ) )/max( 1, norm( evm ) );
         if err > max( 1, n )*tol,
            if Details == 1,  disp( 'Failed 4' ),  end
            failures = failures + 1;
         else
            max_err = max( max_err, err );
         end
         %
         Vx = [ Q Zs(:,1:ns) ];  tolr = tolm*n*max( norm( Vx ), 1 );
         err = max( [ abs( neig - ns ) abs( ns - rank( Vx, tolr ) ) ] );
         if err > 0,  
            if Details == 1,  disp( [ 'Failed 5, orthm = ', num2str( orthm ) ] ),  end
            failures = failures + 1;
         end
      end
   end
end

if max_err < tol,
   disp( [ 'ZGHFDF :    passed  -- maximum relative error max_err = ', num2str( max_err  ) ] )
   disp( [ '            Number of problems solved using randn     = ', num2str( count    ) ] )
   disp( ' ' )
else
   disp( [ 'ZGHFDF :    failed  -- maximum relative error max_err = ', num2str( max_err ) ] )
   disp( [ '            Number of problems solved using randn     = ', num2str( count   ) ] )
end

if failures > 0,  
   disp( [ 'Number of failed tests                                = ', num2str( failures ) ] )
   disp( ' ' )
end        

