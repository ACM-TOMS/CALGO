% test sequence

% preconditions


%% test toolboxes    
lic = license('test', 'symbolic_toolbox');
if( ~lic );
    warning( 'testsequence:symbolic', 'ERROR: Symbolic Toolbox is not installed or licensed. The t-packages may not run.' ); end;

lic = license('test', 'signal_toolbox');
if( ~lic );
    warning( 'testsequence:signal', 'ERROR: Signal Processing Toolbox is not installed or licensed. The t-packages may not run.' ); end;

lic = license( 'test', 'Distrib_Computing_Toolbox' );
pf = isequal( exist('parfor','file'), 2 );
if( ~lic || ~pf );
    warning( 'testsequence:distcomp', 'ERROR: Parallel Toolbox is not installed or licensed. The t-packages may not run.' ); end;



%% test addsequence
EXPECT_EQ( addsequence([1 2],[4;5],[1;1],[1;0]), [4 1 2;5 0 0] );

%% test characteristic
EXPECT_NO_THROW( 'characteristic([1 2 3]'',''amin'',[-1;0;0]);' );
EXPECT_NO_THROW( 'characteristic([1 2 3]'');' );
EXPECT_EQ( characteristic([1 2 3]), [1;1;1] );
EXPECT_EQ( characteristic([1 2 3],'amin',[0]), [0;1;1;1] );
EXPECT_EQ( characteristic([1 2 3],'amin',[0],'amax',5), [0;1;1;1;0;0] );
[a,amin] = characteristic([1 0; 2 0; 1 -1]','amin',[-1; -1]);
EXPECT_EQ( a, [0 0;0 0;1 1;0 1] );
EXPECT_EQ( amin, [-1;-1] );
EXPECT_EQ( characteristic([1 0; 2 0; 1 -1]','amin',[1; -1],'amax',[2;2]), [1 1 0 0;0 1 0 0] );
D1 = [0 0;2 2; 0 1;1 2;0 2]';
[chi,amin,idx] = characteristic( D1 );
D2 = supp( chi, 2, amin ); 
D2 = D2(:,idx);
EXPECT_EQ( D1, D2 );


%% test constructmu
EXPECT_EQ( constructmu(2,2), [2 1 0;0 1 2] );


%% test diffsequence
EXPECT_NO_THROW( 'diffsequence([1 1 1]'',1);' );
val = diffsequence([1 2 3; 4 5 6],[1 0],'cell');
EXPECT_EQ( val{1}, [1 2 3; 3 3 3; -4 -5 -6] );
EXPECT_NO_THROW( 'diffsequence([1 2 3; 4 5 6],[1 0]);' );
EXPECT_NO_THROW( 'diffsequence([1 2 3; 4 5 6],1);' );
EXPECT_EQ( diffsequence([1 2 3; 4 5 6],1,'equalsize'), {[1 2 3 0;3 3 3 0;-4 -5 -6 0]; [1 1 1 -3;4 1 1 -6;0 0 0 0]} );


%% test equalizeminidx
val = equalizeminidx( {[10 ; 10],[1 1]';[1],[0 0]'} );
EXPECT_EQ( val{1}, [0 0;0 10;0 10] );
EXPECT_EQ( val{2}, [1 0;0 0;0 0] );


%% test mask2symbol/symbol2mask
syms x y z z1 z2 a b;
mask = symbol2mask( mask2symbol([1 2 3; 4 5 6]) );
EXPECT_EQ( mask, [1 2 3; 4 5 6] );
EXPECT_NO_THROW( 'mask2symbol([1]);' );
val = symbol2mask(mask2symbol([]));
EXPECT_TRUE( isempty(val) || isequal(val,0) );
EXPECT_EQ( mask2symbol([1 1],'dim',1), 1 + z1 );
EXPECT_EQ( mask2symbol([1 1],'var','x'), 1 + y );
EXPECT_EQ( mask2symbol([1 1],'var',{'b','a'}), 1 + a );
val = mask2symbol([1 1],'amin',[0;1]);
EXPECT_EQ( val, z2^2+z2 );
EXPECT_EQ( symbol2mask(val,'dim',2,'sym','var','z1'), sym([0 1 1]) );


%% test setidx
EXPECT_NO_THROW( 'setidx( [1 2; 2 3], [1 1]'', [0 0]'' );' );
EXPECT_EQ( setidx( [1 2; 2 3], [1 1]', [0 0]', [2 3]' ), [0 0 0 0;0 1 2 0;0 2 3 0] );

%% test supp
EXPECT_EQ( supp([],2,[0]), zeros(2,0) );
EXPECT_EQ( supp([1 1; 0 1],2,[1;-1]), [1 1 2;-1 0 0] );

%% test symbol2mask %tested above

%% test setupsequence %this file


%% test class: sequence
EXPECT_NO_THROW( 'sequence([1 2 3],[1 2]'');' );
if( INIT('all') );
    EXPECT_EQ( ndims(sequence([1 2 3])), 2 );
    EXPECT_EQ( ndims(sequence([1 2 3]')), 1 );
    EXPECT_WARNING( 'sequence([1 2 3 4],[1])' , 'sequence:dim' ); 
end;

c0 = sequence( [1 2 0 4], [0 0] );
c1 = sequence( [-2 2 0 4], [1 1] );
c2 = sequence( [10 20 30 40], [-5; 0] );
c3 = sequence( [1 2 0 4], [1 0] );
c4 = c3;
c4.c(100) = 0;

EXPECT_EQ( c1.c, [-2 2 0 4] );
EXPECT_EQ( c1.idx, [1;1] );
EXPECT_EQ( ref(c1,0), 0 );
EXPECT_EQ( ref(-c1,1), 2 );
EXPECT_EQ( 2*c1, sequence([-4 4 0 8],[1 1]) );
EXPECT_EQ( c1/2, sequence([-1 1 0 2],[1 1]) );
EXPECT_EQ( c1*2, sequence([-4 4 0 8],[1 1]) );
EXPECT_EQ( 2\c1, sequence([-1 1 0 2],[1 1]) );
EXPECT_EQ( c0+c1, sequence([1 2 0 4 0;0 -2 2 0 4],[0 0]) );
EXPECT_EQ( c0+c1, c1+c0 );
EXPECT_EQ( c0.*c1, sequence([0],[0 0]) );
EXPECT_EQ( c0.*c3, c3.*c0 );
EXPECT_NO_THROW( 'c1./c2;' );
EXPECT_EQ( c3<=c1 ,sequence([1 0 1],[1 2]) );
EXPECT_NO_THROW( 'c3<c1;' ); 
EXPECT_NO_THROW( 'c3>c1;' );
EXPECT_NO_THROW( 'c3>=c1;' );
EXPECT_EQ( c3, c4 );
EXPECT_EQ( c3, c3 );
EXPECT_EQ( sum(c1), 4 );
EXPECT_DOUBLE_EQ( norm(c1,-inf), 0 );
EXPECT_DOUBLE_EQ( norm(c1,0), 3 );
EXPECT_DOUBLE_EQ( norm(c1,1), 8 );
EXPECT_DOUBLE_EQ( norm(c1,2), (4+4+16)^(1/2) );
EXPECT_DOUBLE_EQ( norm(c1,3), (8+8+4*4*4)^(1/3) );
EXPECT_DOUBLE_EQ( norm(c1,inf), 4 );
EXPECT_EQ( supp(c1), [1 1 1;1 2 4] );
EXPECT_EQ( c1.idxmin, [1;1] );
EXPECT_EQ( c1.idxmax, [1;4] );
EXPECT_EQ( nnz(c1), 3 );
EXPECT_EQ( ndims(c1), 2 );
EXPECT_NO_THROW( 'numel(c1);' );
EXPECT_NO_THROW( 'size(c1);' );
[cc1,cc2] = equalize(c1,c2);
EXPECT_NO_THROW( 'cc1+cc2;' );
EXPECT_EQ( shrink(cc1), c1 );


EXPECT_EQ( shrink(sequence([0 0 1 0 1],[1 0])), sequence([1 0 1],[1 2]) );

EXPECT_NO_THROW( 'symbol(c1);' );
EXPECT_EQ( upsample(c1,[2 1; 0 2]), sequence([-2 0 0 0 0 0 0;0 0 2 0 0 0 0;0 0 0 0 0 0 0;0 0 0 0 0 0 4],[3 2]) );
EXPECT_EQ( conv(c1,c2,c1), sequence([40 0 0 -160 -360 0 0 960 480 640],[-3 2]) );

EXPECT_EQ( characteristic(c1), sequence([1 1 0 1],[1 1]) );
val = diffsequence(c1,[1 1]);
EXPECT_EQ( val, sequence([-2 4 -2 4 -4;2 -4 2 -4 4],[1 1]) );
