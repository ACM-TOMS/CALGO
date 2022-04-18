%#ok<*SAGROW>
%#ok<*NOPTS>
 %#ok<*NASGU>

%% Example 4.3
S = getS( 'a',sym(1/12)*[3 3 4 3 3 4 3 3 4 3 3].', 'M',-3, 'D',[-2 -1 0] )
[T,Om] = transitionmatrix( S )
V = constructVt( Om, 0 ) 
TV = restrictmatrix( T, V )
vdisp( -12.*TV )

tjsr( TV ); %with balancing
tjsr( TV, 'nobalancing' );

%% Example 4.4
E1 = [2 1;-1 2];
E2 = [2 0; 2 1];
tjsr( {E1,E2}, 'ordering',{[1 2]' [2 0]'}, 'smpflag',[0 1], ...
   'balancingvector',[1 0.95], 'invariantsubspace','none' )

%% Table 1 / Table 2
%  The results are obtained with scripts similar to
dim = 10;
T = [];
i = 0;
while( i<10);
    M = tgallery( 'rand_gauss', dim, 2, 'norm' );
    [r,info] = tjsr( M );
    if( numel(r)==1 && info.info.errorcode<0 );
        i = i+1; 
        T(end+1) = info.counter.treetime; end;  end; 
mean(T)

%% Table 3
%  The results are obtained with scripts similar to
dim = 2;
J = 4;
i = 0;
[GRIP,MODGRIP] = deal( [] );
while( i<3 );
    M = tgallery( 'rand_gauss', dim, J, 'rho' );
    [r,info] = tjsr( M );
    if( numel(r)~=1 || info.info.errorcode>=0 );
        fprintf( ' Modified Invariant polytope algorithm did not find s.m.p..\n' );
        continue; end;
    i=i+1;
    
    [~,~,val] = findsmpold( M, 'gripenberg' );
    if( abs(val.jsrbound(1)-r)<1e-12 )
        GRIP(end+1) = val.time;
    else
        GRIP(end+1) = inf; end;
    
    [~,~,val] = findsmpold( M );
    if( abs(val.jsrbound(1)-r)<1e-12 )
        MODGRIP(end+1) = val.time;
    else
        MODGRIP(end+1) = false; end; end;

%% Example 5.1 / Table 4 / Table 5
%  The results are obtained with scripts similar to
X = tgallery('mejstrik_119');
cmodgrip = findsmpold( X, 'maxsmpdepth',120 )
cgrip = findsmpold( X, 'gripenberg', 'delta',1, 'maxsmpdepth',120, 'v',2 )
cgen = findsmpold( X, 'genetic' )

C15 = tgallery('mejstrik_Cn',15)
cmodgrip = findsmpold( C15, 'maxsmpdepth',20 )

%% Table 6 / Table 7
%  The results are obtained with scripts similar to
D = {[1 1],[1 -1],[-1 1],[-1 -1]};
C = codecapacity( D );
tjsr( C, 'v',2, 'maxsmpdepth',2 )
