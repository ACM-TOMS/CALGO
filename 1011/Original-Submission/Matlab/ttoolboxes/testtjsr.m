% test sequence

%#ok<*CTPCT>
%#ok<*ASGLU>
%#ok<*NASGU>

if( exist('blf','file')~=0 ); 
    subdivisioninstalled = 1; 
else; 
    subdivisioninstalled = 0;  end;


try
    evalc( 'gurobi_setup' ); % run evalc to supress output
    model.A = sparse( [1 1 0; 0 1 1] );
    model.obj = [1 2 3];
    model.modelsense = 'Max';
    model.rhs = [1 1];
    model.sense = [ '<' '<'];
    params.outputflag=0;

    gurobiresult = gurobi( model, params );
    gurobiinstalled = true;
    if( ~isequal(gurobiresult.x,[ 1 0 1].') ); 
        fprintf( 'ERROR: Gurobi Broken\n ==========================\nThe tjsr-package will NOT work and may NOT emit error-messages.\n' ); end; 
catch
    gurobiinstalled = false; end;


inst = ~isempty( ver('distcomp') );
lic = license( 'test', 'Distrib_Computing_Toolbox' );
file = isequal( exist('parfor','file'), 2 );
if( ~lic || ~file || ~inst );
    parinstalled = false;
else
    parinstalled = true;
end;

verbose = 0;



% preconditions    

%% test external functions / libraries

%% test version
EXPECT_FALSE( verLessThan('matlab','9.1') ); 

%% test binarymatrix
val = binarymatrix( 2, 2, 2^5+2^3+2^2+2^1 );
EXPECT_EQ( val, {[0 1; 0 0],[1 1; 1 0]} );
if( INIT('all') )
    EXPECT_PRED( @isempty, binarymatrix( 1, 1, 2 ) );
    
    EXPECT_PRED( @isempty, binarymatrix(3, 1, 1+16+32+128+256,1) );
    EXPECT_EQ( binarymatrix(3, 1, 1+16+32+128+256,0), {[1 1 0;1 1 0;0 0 1]} );
    EXPECT_NPRED( @isempty, binarymatrix(3, 1, 1+2+8+16+256,0) );
    EXPECT_PRED( @isempty, binarymatrix(2, 1, 4,1) );
    EXPECT_NPRED( @isempty, binarymatrix(2, 1, 2,1) );
end

%% test blockjsr
EXPECT_EQ( blockjsr([1 0],[.5 2]), [.5 2] );
EXPECT_EQ( blockjsr(1), 1 );
EXPECT_EQ( blockjsr(0,1), 1 );
EXPECT_EQ( blockjsr([0 1],1), 1 );
if( INIT('all') );
    EXPECT_ERROR( 'blockjsr([1 2 3]);' );
end

%% test chooseval
EXPECT_EQ( [1 0 0], chooseval([inf 0 nan],1) );
if( INIT('all') );
    EXPECT_PRED( @islogical, chooseval(randn(1,100),[10 50]) );
    EXPECT_PRED( @islogical, chooseval(randn(1,100),[10 50]) );
    EXPECT_EQ( chooseval([0 0 0 0 3 3]'), [0 0 0 0 1 1]' );
    EXPECT_EQ( chooseval([1 1 1 4 4 4]',[0 0]), [0 0 0 0 0 0]');
end

%% test codecapacity
EXPECT_NO_THROW( 'codecapacity( {[1 1],[-1 -1]}, ''v'',verbose );' );
if( INIT('all') );
    EXPECT_NO_THROW( 'codecapacity( {[1 1 0],[1 -1 0]}, ''v'',verbose, ''plot'',1 ); title( ''codecapacity'' );' );
    C = codecapacity( {[1 1],[-1 -1]}, 'v',verbose );
    EXPECT_EQ( C, {[0 1; 1 1],[1 1; 1 0]} );
end


%% test computepolytopenorm
EXPECT_NO_THROW( 'val = computepolytopenorm( [0 3]'', [], 1, 0, 0, 0 );' );
EXPECT_EQ( val, inf );
if( INIT('all') );
    if( gurobiinstalled && parinstalled )
        EXPECT_PRED( @isempty, computepolytopenorm( [], [1 2; 2 1], 1, 1, 0, 0, 'g' ) );
        EXPECT_NEAR( computepolytopenorm( [0 3; 2 2]', [1 2; 2 1]', 1, 0, 0, 0, 'g' ), [3 4/3], 5e-12 );
        EXPECT_DOUBLE_EQ( computepolytopenorm([0 3; 2 2]', [1 2; 2 1]', 0, 0, 0, 0, 'g' ), [3/2 4/3] );
    end;
    EXPECT_PRED( @isempty, computepolytopenorm( [], [1 2; 2 1], 1, 1, 0, 0, 'm' ));
    EXPECT_DOUBLE_EQ( computepolytopenorm( [0 3; 2 2]', [1 2; 2 1]', 1, 0, 0, 0, 'm' ), [3 4/3] );
    EXPECT_DOUBLE_EQ( computepolytopenorm( [0 3; 2 2]', [1 2; 2 1]', 0, 0, 0, 0, 'm' ), [3/2 4/3] );
end

%% test daubechiesmatrix
if( subdivisioninstalled );
    EXPECT_NO_THROW( 'daubechiesmatrix(3,''gugl'',''v'',verbose,''double'');' );
    D21 = daubechiesmatrix( 2, 'gugl', 'v',verbose ); 
    val1 = max(cellfun(@(x) max(abs(eig(double(x)))), D21) );
    D22 = daubechiesmatrix( 2,'jung','v',verbose); 
    val2 = max( cellfun(@(x) max(abs(eig(double(x)))), D22) );
    EXPECT_DOUBLE_EQ( val1, val2 );
end

%% test estimatejsr
TT = {[1 2; 0 1],[1 2; 1 0]};
EXPECT_NO_THROW( 'val =  estimatejsr( TT );' );
EXPECT_NEAR( val(1), 2, 1e-9 );
EXPECT_LE( val(2), 2.5 );


%% test estimatepolytopenorm
VV = [ -6    11    -8    -7     8     3     3     3    -4     5
        7   -17     3    -1   -15    10    -6    21    -3    13
       -3     8     1    13     2   -19     3    -9   -15    -8];
pp = [-10    -3     7     9    -6    25    -7    12    -3    -7
        4    -3    -5    13     3    14     3     0    -6     8
       12    -2    -8    18   -16   -13    -8     6   -13     0];

epsilon = 1e-10;
EXPECT_NO_THROW( 'estimatepolytopenorm(pp, VV, [], 1, epsilon);' );
if( INIT('all') )
    [normval10,bd10] = estimatepolytopenorm( abs(pp), abs(VV), [], 10, epsilon );
    [normval,bd] = estimatepolytopenorm( abs(pp), abs(VV), [], 0, epsilon );
    idxin = normval<bd(1); 
    idxout = normval>bd(2);
    EXPECT_FALSE( any([normval10(idxin)>1  normval10(idxout)<1]) );

    [normval11,bd11] = estimatepolytopenorm( pp, VV, [], 11, epsilon );                        
    [normval,bd] = estimatepolytopenorm( pp, VV, [], 1, epsilon );
    EXPECT_FALSE( any([normval11(normval<bd(1))>1  normval11(normval>bd(2))<1]) );

    %[normval,bd] = estimatepolytopenorm( pp, VV, [], 55, epsilon ); % XX experimental
    %EXPECT_FALSE( any([normval11(normval<bd(1))>1  normval11(normval>bd(2))<1]) );

    [normval,bd] = estimatepolytopenorm( pp, VV, [], 6, epsilon ); 
    EXPECT_FALSE( any([normval11(normval<bd(1))>1  normval11(normval>bd(2))<1]) );


    [normval,bd] = estimatepolytopenorm( pp, VV, [], 7, epsilon ); 
    EXPECT_FALSE( any([normval11(normval<bd(1))>1  normval11(normval>bd(2))<1]) );

    VV2 = [-6    11    -8    -7     8     3     3     3    -4     5
            7   -17     3    -1   -15    10    -6    21    -3    13];
    pp2 = [-10   -3     7     9    -6    25    -7    12    -3    -7
             4   -3    -5    13     3    14     3     0    -6     8];        

    EXPECT_NO_THROW( '[normval11,bd11] = estimatepolytopenorm(pp2, VV2, [], 11, epsilon);' );
    [normval,bd] = estimatepolytopenorm( pp2, VV2, [], 5, epsilon ); 
    EXPECT_FALSE( any([normval11(normval11<bd(1))>1  normval11(normval>bd(2))<1]) );

end


%% test extravertex
V = [1 2 3;2 3 4]';
[Ve,VVe] = extravertex( V );
EXPECT_EQ( [V,Ve],VVe );
EXPECT_GE( min(svd(VVe)), .01 );
if( INIT('all') )
    [Ve,VVe] = extravertex( randn(10), 't', 2000 );
    EXPECT_EQ( size(Ve,2),10 );

    [Ve,VVe] = extravertex( randn(10), 't', 0 );
    EXPECT_EQ( size(VVe,2),10 );

    [Ve,VVe] = extravertex( V, 'pt', 0 );
    EXPECT_PRED( @isempty, Ve  );

    EXPECT_EQ( eye(3), extravertex( zeros(3,0), 'pt', 0 ) );
    EXPECT_EQ( eye(3), extravertex( zeros(3,0), 'pt', 1 ) );

    V2 = [1 1 0]';
    [Ve,VVe] = extravertex( V2, 'pt', 0 );
    EXPECT_EQ( Ve/max(Ve), [0;0;1] );

    [Ve,VVe] = extravertex( V2, 'pt', 0, 'c' );
    EXPECT_PRED( @iscell, Ve );
    EXPECT_PRED( @iscell, VVe );

    M = {[1 2;1 2], [2 4;2 4]};
    Mb = {[2 4;2 4], [4 8;4 8]};
    oo = {[1 2 2 1 1 2 1]'};
    v0 = {[1 0]'};
    Ve1 = extravertex( M, oo, v0 );
    Ve1a = extravertex( M, oo, v0 ,'s', 2 );
    Ve1b = extravertex( Mb, oo, v0 );
    EXPECT_EQ( Ve1a, Ve1b );
    EXPECT_EQ( size(Ve1,2), 1 );

    M = {[1 2;1 2], [2 4;2 4]};
    oo = [1 2 2 1 1 2 1]';
    v0 = [1 0]';
    Ve2 = extravertex( M, oo, v0 );
    EXPECT_EQ( size(Ve2,2), 1 );
    EXPECT_EQ( Ve1, Ve2 );

    VV = [10000 0 0;10000 0 1;10000 1 0;10000 1 1]';
    Ve = extravertex( VV );
    EXPECT_EQ( size(Ve,2), 2 );
end



%% test findsmp
AA = {[1 -1; 3 -2], [1 3; -1 -1]};
EXPECT_NO_THROW( 'findsmp( AA, ''v'',0 );' );
if( INIT('all') )
    fprintf( '\nfindsmp: ........... (This may take long)\n         ' );
    AA = {[1 -1; 3 -2], [1 3; -1 -1]};
    fprintf( '|' ); 
    findsmp( AA, 'bruteforce', 'minsmpdepth',3, 'v',verbose, 'b',[0 inf] );
    findsmp( AA, 'nlbf', 'minsmpdepth',7, 'maxsmpdepth',8, 'v',verbose );
    findsmp( AA, 'gripenberg', 'maxsmpdepth',5, 'v',verbose );
    findsmp( AA, 'genetic', 'v',-1 );

    BB = tgallery( 'mejstrik_Cn', 30 );
    [c4,nc,info] = findsmp( BB, 'maxsmpdepth',inf, 'v',verbose, 'bound', [0 1.03] );
    EXPECT_LE( length(c4{1}), 25  );
    EXPECT_LT( info.jsrbound(1), 1.0338951135135 );

    [c4,nc,info] = findsmp( BB, 'maxsmpdepth',inf, 'v',verbose, 'bound', [1.01 1.2] );
    EXPECT_LE( length(c4{1}), 20 );
    EXPECT_LT( info.jsrbound(1), 1.019 ); 
    EXPECT_GT( info.jsrbound(2), 1.197 );

    %
    AA = {[1 -1; 3 -2], [1 3; -1 -1]};
    fprintf( '|' ); 
    EXPECT_EQ( findsmp( AA, 'modgrip', 'maxsmpdepth',5, 'v',verbose ), {[1;2]} );
    EXPECT_EQ( findsmp( AA, 'modgrip', 'maxsmpdepth',5, 'v',verbose, 'vpa' ), {[1;2]} );
    EXPECT_EQ( findsmp( AA, 'lowgrip', 'maxsmpdepth',5, 'v',verbose ), {[1;2]});
    EXPECT_EQ( findsmp( AA, 'lowgrip', 'maxsmpdepth',5, 'v',verbose, 'vpa' ), {[1;2]} );
    EXPECT_EQ( findsmp( AA, 'randgrip', 'maxsmpdepth',5, 'v',verbose ), {[1;2]} );
    EXPECT_EQ( findsmp( AA, 'randgrip', 'maxsmpdepth',5, 'v',verbose, 'vpa' ), {[1;2]} );
    %
    AA = tgallery( 'rand_gauss', 100, 2, 'sparse',.8, 'rho' );
    fprintf( '|' ); 
    findsmp( AA, 'maxsmpdepth',10, 'v',verbose );
    %
    AA = tgallery('rand_gauss',30,2,'seed',102,'rho');    
    fprintf( '|' ); 
    findsmp( AA, 'maxsmpdepth',20, 'v',verbose, 'nearlycandidate',.99 );
    %
    AA = tgallery( 'rand_gauss', 10, 100, 'seed',102, 'rho' );
    fprintf( '|' ); 
    findsmp( AA, 'modgrip', 'maxsmpdepth',20, 'v',verbose, 'N',2 );
    findsmp( AA, 'randgrip', 'maxsmpdepth',20, 'v',verbose, 'N',2 );
    findsmp( AA, 'highgrip', 'maxsmpdepth',20, 'v',verbose, 'N',2 );
    %
    AA = tgallery( 'rand_gauss', 10, 3, 'seed',102, 'rho' );    
    findsmp( AA, 'randgrip', 'maxsmpdepth',5, 'v',verbose, 'N',inf );
    %
    fprintf( '|' );
    EXPECT_EQ( findsmp(daubechiesmatrix( 7, 'v',0, 'd' ), 'v',verbose),{[1],[2]} );
    fprintf( '|' );
    EXPECT_EQ( findsmp(daubechiesmatrix( 2, 'v',0, 'd' ), 'v',verbose),{[1]} );
    fprintf( '|' );
    EXPECT_EQ( findsmp(daubechiesmatrix( 10, 'v',0, 'd' ), 'v',verbose),{[1;1;2;2]} );
    fprintf( '|' );
    EXPECT_EQ( findsmp(tgallery('grip_p52'), 'v',0 ),{[1;1;1;1;1;1;1;1;1;1;1;1;2]} );
    fprintf( '|' );
    EXPECT_EQ( findsmp(tgallery('mejstrik_Cn',10),'v',0),{[1;1;1;1;1;1;1;1;1;1;2]} );
    fprintf( '|' );
    TT = {[1 2; 0 1],[1 2; 1 0]};
    [c,nc,i] = findsmp( TT, 'nc',.4, 'v',0 );
    EXPECT_FALSE( ~isequal(c,{2}) || ~isequal(nc,{1}) || min(i.jsrbound)<2-1e-10 || max(i.jsrbound)>2.1 );
    
    fprintf( '\n' ); 
end


%% test intersectinterval
EXPECT_EQ( intersectinterval([4 5],[1 4],[-inf 4]), [4 4] );
if( INIT('all') )
    EXPECT_EQ( intersectinterval(), [-inf inf] );
    EXPECT_EQ( intersectinterval([0 0]), [0 0] );
    EXPECT_EQ( intersectinterval([0 0],[1 1]), [] );
    EXPECT_EQ( intersectinterval([-1 0],[1 2],'failsafe','minmax'), [-1 2] );
    EXPECT_EQ( intersectinterval([-1 0],[1 2],'failsafe','maxmin'), [0 1] );
    EXPECT_EQ( intersectinterval([1 2],2), 2 );
    EXPECT_EQ( intersectinterval([1 2],1) ,1 );
    EXPECT_EQ( intersectinterval([0 5],1,[1 3]), 1 );
    EXPECT_EQ( intersectinterval([0 5],6), [] );
    EXPECT_EQ( intersectinterval([0 5],[6 7]), [] );
    EXPECT_ERROR( 'intersectinterval([-1 0],[1 2],''failsafe'',''wrongoption'');' );
end

%% test invariantsubspace
[M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'v',verbose );
EXPECT_EQ( size(M,2),3);
if( INIT('all') )
    [M,B] = invariantsubspace( {sym([1 1 1; 1 0 1; 0 1 0]), sym([1 -1 0; 0 2 0; 1 1 2])}, 'v',verbose );
    EXPECT_EQ( numel(M), 3 );
    [M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'sym', 'v',verbose );
    EXPECT_EQ( numel(M), 3 );
    [M,B] = invariantsubspace( {sym([1 1 1; 1 0 1; 0 1 0]), sym([1 -1 0; 0 2 0; 1 1 2])}, 'double', 'v',verbose );
    EXPECT_EQ( numel(M), 3 );
    [M,B] = invariantsubspace( {sym([1 1 1; 1 0 1; 0 1 0]), sym([1 -1 0; 0 2 0; 1 1 2])}, 'vpa', 'v',verbose );
    EXPECT_EQ( numel(M), 3 );
    [M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'vpa', 'v',verbose );
    EXPECT_EQ( numel(M), 3 );

    [M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'vpa', 'v',verbose,'none' );
    EXPECT_EQ( B,eye(3) );
    [M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'trans', 'd', 'v',verbose-1 );
    EXPECT_EQ( numel(M), 2 );
    [M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'perm', 'd' );
    EXPECT_EQ( numel(M), 1 );
    [M,B] = invariantsubspace( {[1 1 1; 1 0 1; 0 1 0], [1 -1 0; 0 2 0; 1 1 2]}, 'perm', 'd', 'sym' );
    EXPECT_EQ( numel(M), 1 );
end

%% leadingeigenvector
MM = [1 2;1 1];
[M,v0,v0s,mult] = leadingeigenvector( MM, 'v',verbose );
EXPECT_EQ( M{1}, MM );
EXPECT_DOUBLE_EQ( v0{1}, 1/sqrt(3)*[sqrt(2);1] );
EXPECT_DOUBLE_EQ( v0s{1}, [1;sqrt(2)]/dot(v0{1},[1;sqrt(2)]) );
EXPECT_EQ( mult, 1 );
    
if( INIT('all') )
    [M,v0,v0s,mult] = leadingeigenvector( [1 2;-2 1], 'v',verbose, 'cp',0 );
    EXPECT_EQ( numel(M),2);

    [M,v0,v0s,mult] = leadingeigenvector( {[1 1;-1 1], [1 0;1 2]}, 'v',verbose, 'cp',0 );
    v0 = cell2mat( v0 );
    v0 = normalizematrix( v0, 'colnorm',inf ); 
    v0 = sortrows( v0.' ).';
    EXPECT_EQ( size(v0), [2 3] );
    EXPECT_DOUBLE_EQ( v0, [0 1;1 -1i;1 1i].' );

    [M,v0,v0s,mult] = leadingeigenvector( [1 2 0;-2 1 0; 0 0 sqrt( 5 )], 'v',verbose, 'cp',1 );
    EXPECT_FALSE( ~isequal(M{1},[1 2 0;-2 1 0; 0 0 sqrt( 5 )]) || ~isequal(v0{1},[0;0;1]) || ~isequal(v0s{1},[0;0;1]) || ~isequal(mult,1) );  
        

    [M,v0,v0s,mult] = leadingeigenvector( {[0 1;-1 0], [1 0;0 1]}, 'v',verbose, 'cp',1 );
    v0 = cell2mat( v0 );
    v0 = normalizematrix( v0, 'colnorm',inf ); 
    v0 = sortrows( v0.' ).';
    EXPECT_FALSE( ~isequal(size(v0),[2,4]) || norm(v0-[0 1;1 0;1 -1i;1 1i].')>1e-12 );  

    [M,v0,v0s,mult] = leadingeigenvector( [1 2 0;-2 1 0; 0 0 sqrt( 5 )], 'v',verbose, 'cp',2 );
    EXPECT_FALSE( numel(v0)~=1 || norm(v0{1}-[0;0;1])>1e-12 );

    [M,v0,v0s,mult] = leadingeigenvector( [0 1;-1 0], 'v',verbose, 'cp',2 );
    v0 = cell2mat( v0 );
    v0 = normalizematrix( v0, 'colnorm',inf ); 
    v0 = sortrows( v0.' ).';
    EXPECT_FALSE( ~isequal(size(v0),[2 2]) || norm(v0-[1 -1i;1 1i].')>1e-12 );

    [M,v0,v0s,mult] = leadingeigenvector( {[0 1;-1 0],[1 0;0 1]}, 'v',verbose, 'cp',2 );
    v0 = cell2mat( v0 );
    v0 = normalizematrix( v0, 'colnorm',inf ); 
    v0 = sortrows( v0.' ).';        
    EXPECT_EQ( size(v0), [2 2] );
    EXPECT_NEAR( v0, [0 1;1 0], 1e-12 );

    [M,v0,v0s,mult] = leadingeigenvector( [1 0;0 1], 'v',verbose, 'cp',3 ); %XX This option is not correctly implemented yet
    v0 = cell2mat( v0 );
    v0 = normalizematrix( v0, 'colnorm',inf ); 
    v0 = sortrows( v0.' ).';
    EXPECT_EQ( size(v0), [2 2] );
    EXPECT_NEAR( v0, [0 1;1 0], 1e-12 );


    [M,v0,v0s,mult] = leadingeigenvector( [1 2 0;-2 1 0; 0 0 sqrt( 5 )], 'v',verbose, 'cp',4 );
    EXPECT_EQ( numel(v0), 1 );
    EXPECT_NEAR( v0{1}, [0;0;1], 1e-12 );

    [M,v0,v0s,mult] = leadingeigenvector( [0 1;-1 0], 'v',verbose, 'cp',4 );
    EXPECT_EQ( numel(v0), 0 );

    TT = {[1 2 ; 1 0],[0 0; 1 0]}; 
    oo = {[1 2]',[1]'}; %#ok<NBRAK>
    [oret,v0,v0s,mult,oorig] = leadingeigenvector( TT, oo, 'v',verbose, 'c',0, 'r',0 );
    v0 = cell2mat( v0 );
    v0 = normalizematrix( v0, 'colnorm',inf ); 
    v0 = sortrows( v0.' ).';
    EXPECT_EQ( oret, {[1;2],1} );
    EXPECT_NEAR( v0, [0 1;1 .5]', 1e-12 );

    [oret,v0,v0s,mult,oorig] = leadingeigenvector( TT, oo, 'v',verbose, 'c',1, 'r',50 );
    %EXPECT_EQ( oret, {[1;2],[2;1],1,[1;1]} );

    %[oret,v0,v0s,mult,oorig] = leadingeigenvector( TT, oo, 'v',verbose );
    %EXPECT_EQ( oret, {[1;2],[2;1],1} );
end;



%% test makeorderinggraph
oo = [1 1 1 2 1 2 2 2 1 1;0 2 1 0 2 1 2 1 1 2;0 0 0 0 1 0 0 1 2 1]; 
EXPECT_NO_THROW( 'makeorderinggraph( oo, ''labeledge'',0, ''v'',verbose );' );


%% testpartitionatepolytope
EXPECT_NO_THROW( 'partitionatepolytope(randn(2,100),2,0);' );

%% test preprocessmatrix
EXPECT_NO_THROW( 'preprocessmatrix({[-1 2; 2 3],[-1 2; 2 3]},''v'',verbose);' );
EXPECT_EQ( preprocessmatrix({-1}, 'v',verbose), {1} );
if( INIT('all') )
    EXPECT_NO_THROW( 'preprocessmatrix({[-1 2; 2 3],[-1 2; 2 3]},''v'',verbose,''inverse'',1,''addinverse'',1,''transpose'',1,''addtranspose'',1,''makepositive'',0,''timestep'',.9,''perturbate'',.001,''removezero'',0,''removeduplicate'',0,''basechange'',''random'',''exponential'',1);');
    EXPECT_NO_THROW( 'preprocessmatrix({[-1 2; 2 3],[-1 2; 2 3]},''v'',verbose,''basechange'',1);');
    EXPECT_NO_THROW( 'preprocessmatrix({[-1 2; 2 3],[-1 2; 2 3]},''v'',verbose,''basechange'',[1 1;1 -1]);');
end

%% test reducelength
val = reducelength( {[2 1 2 2 1 2],[3 1]} );
EXPECT_EQ( val, {[1 2 2],[1 3]} );
if( INIT('all') )
    EXPECT_EQ( reducelength({[]},0),{[]} );
    EXPECT_EQ( reducelength({[2 1 2 2 1 2],[3 1]},0),{[1 2 2],[1 3]} );
    EXPECT_EQ( reducelength([2 1 2 2 1 2],1),[2 1 2] );
    EXPECT_EQ( reducelength([2 1 2 2 1 2],2),[1 2 2 1 2 2] );
end

%% test removecombination
val = removecombination( {[1 1 2],[1]} );
EXPECT_FALSE( size(val,2)~=2 || ~searchincellarray(2,val,1) || ~searchincellarray(1,val,1) );
if( INIT('all') )
    EXPECT_EQ( removecombination({[1 1 2],[1]},'add'), {[1 1 2],[1],[1 2],[2]} );
end

%% test tbuildproduct_fast
EXPECT_EQ( tbuildproduct_fast({2 4 5},[1 0 3]), 10 );
if( INIT('all') )
    vv = [1;1];
    val1 = tbuildproduct_fast( {[2 1; 1 2],[2 0;0 1]}, [1 0 2] );
    val2 = tbuildproduct_fast( {[2 1; 1 2],[2 0;0 1]}, [1 0 2], vv );
    EXPECT_DOUBLE_EQ( val1*vv, val2 );
end

%% test testtjsr

%% test tgallery
T = tgallery( 'rand', 5, 2 );
if( INIT('all') )
    EXPECT_PRED( @(x) identifymatrix(x,0,'columnstochastic'), tgallery('rand_stochastic',5,2) );
    
    EXPECT_PRED( @(x) identifymatrix(x,0,'doublestochastic'), tgallery('rand_doublestochastic',5,2) );
    EXPECT_NPRED( @isempty, tgallery('rand_stochastic_neg',5,2)  );
    EXPECT_NPRED( @isempty, tgallery('rand_doublestochastic_neg',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,0,'pm1'), tgallery('rand_pm1',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,0,'bool'), tgallery('rand_bool',5,2) );
    EXPECT_NPRED( @isempty, tgallery('rand_gauss',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,1,'sparsity'), tgallery('rand_gauss',5,2,'sparse',.5) );
    EXPECT_PRED(  @(x) identifymatrix(x,1,'int'), tgallery('rand_gauss',5,2,'int',1) );
    EXPECT_PRED( @(x) identifymatrix(x,0,'pos'), tgallery('rand',5,2) );
    EXPECT_PRED(  @(x) identifymatrix(x,0,'pos'), tgallery('rand_equal',5,2));
    EXPECT_NPRED( @isempty, tgallery('rand_neg',5,2) );
    EXPECT_NPRED( @isempty, tgallery('rand_equal_neg',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,1,'unitary'), tgallery('rand_unitary',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,1,'zero'), tgallery('rand_zero',5,2) );
    EXPECT_NPRED( @isempty, tgallery('rand_colu_1',5,2) );
    EXPECT_NPRED( @isempty, tgallery('rand_colu_0',5,2));
    EXPECT_NPRED( @isempty, tgallery('rand_corr_1',5,2));
    EXPECT_NPRED( @isempty, tgallery('rand_corr_0',5,2));
    EXPECT_PRED( @(x) identifymatrix(x,0,'hessu'), tgallery('rand_hess',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,1,'unitary'), tgallery('orthog_1',5,2) );
    EXPECT_PRED( @(x) identifymatrix(x,1,'unitary'), tgallery('orthog_2',5,2));
    EXPECT_PRED( @(x) identifymatrix(x,1,'unitary'), tgallery('orthog_3',5,2));
    EXPECT_NPRED( @isempty, tgallery('cex') );
    EXPECT_NPRED( @isempty, tgallery('binary',5,2,10) );
    EXPECT_PRED( @isempty, tgallery('binary2',5,2,10) );
    EXPECT_NPRED( @isempty, tgallery('mejstrik_119',5,2,10) );
    EXPECT_NPRED( @isempty, tgallery('mejstrik_Cn',5,2,10) );
    EXPECT_NO_THROW( 'tgallery(''rand_gauss'',5,2,''seed'',rng,''sparse'',.5,''bool'');' );
end

%% test tjsr
EXPECT_NO_THROW( 'tjsr( {[1 1;1 0],[1 -1; 1 1]}, ''v'',0 );' );
if( INIT('all') )
    %XX Add tests which also test the return value

    fprintf( '\ntjsr:................... (This will take long)\n     ' ); 
    %
    fprintf( '|' ); 
    tjsr( tgallery('rand_pm1',3,2,'seed',100), 'ordering',{[1 2 2 2 2 2 2 2]'}, 'plot','info_rho', 'v',verbose-1 );
    EXPECT_EQ( tjsr({[1],[2]},'v',verbose-1), 2 );
    EXPECT_NEAR( tjsr({[1 -1; 0 1],[0 0; 1 0]}, 'v',verbose-1), 1.31950791077, 3e-12 );
    %
    fprintf( '|' );
    AA = {[1 1 0; 1 0 1; 1 0 0],[1 0 1; 1 1 0; 0 1 0]};
    [JSR,info] = tjsr( AA, 'plot','polytope', 'proof',1, 'v',-1 );
    EXPECT_NEAR( JSR, 1.83928675521, 5e-12 );
    EXPECT_EQ( info.info.errorcode, -10 );
    EXPECT_EQ( numel(JSR), 1 );
    %
    % stupid input tests
    fprintf( '|' );
    clf;
    B = {[1 1 0; -1 0 -1; 1 0 0],[1 0 -1; 1 -1 0; 0 1 0]};
    EXPECT_NEAR( tjsr( B, 'v',verbose-1 ), 1.39406936116, 2e-12);
    EXPECT_PRED( @isnan, tjsr({nan},'v',verbose-1) );
    EXPECT_PRED( @isinf, tjsr({inf},'v',verbose-1) );
    EXPECT_EQ( tjsr({[1],[2]},'v',verbose-1), 2 );
    EXPECT_EQ( tjsr({[1],[-1]},'v',verbose-1), 1 );
    EXPECT_EQ( tjsr({[1 2; 1 2]},'v',verbose-1), 3 );
    %
    AA = {[1 1;0 1],[1 2;0 2]}; 
    EXPECT_EQ( tjsr(AA,'v',verbose-1), 2 );
    %
    AA = {[1 1 1;0 1 -1;1 -1 -1], [0 0 1;-1 -1 0;1 1 -1]};
    fprintf( '|' ); JSR = tjsr( AA, 'ordering',{[2]'}, 'v', verbose-1, 'maxiteration',10 ); 
    EXPECT_NEAR( JSR, 2, 10e-12 );
    EXPECT_EQ( numel(JSR), 1 );
    fprintf( '|' ); JSR = tjsr( AA, 'extravertex', {[.1 0 0]'},'v', verbose-1,'invariantsubspace','none', 'maxiteration',10 );  
    EXPECT_NEAR( JSR, 2, 10e-12 );
    EXPECT_EQ( numel(JSR), 1 );
    fprintf( '|' ); JSR = tjsr( AA, 'ordering', {[1 2 2 2 2 2 2 2 2]',[]'}, 'smpflag',[0 2], 'v0', {[0.058585928823279 -0.687551547968790 0.723768303968635]',[.1 0 0]'}, 'v',verbose-1, 'invariantsubspace','none', 'maxiteration',10 ); 
    EXPECT_NEAR( JSR, 2, 10e-12 );
    EXPECT_EQ( numel(JSR), 1 );
    %
    D7 = daubechiesmatrix( 7, 'v',0 );
    fprintf( '|' ); JSR = tjsr( D7, 'v',verbose-1, 'maxiteration',10 ); 
    EXPECT_NEAR( JSR, 0.181695060514689, 10e-12 );
    EXPECT_EQ( numel(JSR), 1 );
    fprintf( '|' ); JSR = tjsr( D7, 'balancingvector',[1 1.02 .01 .01 .01 .01 .01], 'v',verbose-1, 'maxiteration',10, 'invariantsubspace','none' ); 
    EXPECT_NEAR( JSR, 0.181695060514703, 10e-12 );
    EXPECT_EQ( numel(JSR), 1 );
    fprintf( '|' ); JSR = tjsr( D7, 'nobalancing',1, 'v', verbose-1, 'maxiteration',100, 'epspolytope',-1e-8 ); 
    EXPECT_GE( numel(JSR), 2 );
    fprintf( '|' ); JSR = tjsr( D7, 'nobalancing',1, 'delta',.99999, 'v',verbose-1, 'maxiteration',13 ); 
    EXPECT_NEAR( JSR(1), 0.181695060514703, 10e-4 );
    EXPECT_NEAR( JSR(2), 0.181695060514703, 10e-4 );
    EXPECT_EQ( numel(JSR), 2 );
    fprintf( '|' ); JSR = tjsr( D7, 'alwaysout', 'maxiteration',3, 'v',verbose-1, 'fastnorm',0 );
    %
    % plot tests and other stuff
    T = tgallery( 'rand_gauss', 20, 2, 'rho', 'seed',100, 'v',verbose-1 );
    fprintf( '|' ); tjsr( T, 'naturalselectiontype',2, 'plot','norm', 'v',verbose-1, 'maxiteration',4 );
    fprintf( '|' ); figure; tjsr( T, 'naturalselectiontype',3, 'plot','L', 'v',verbose-1, 'maxiteration',5 );
    fprintf( '|' ); tjsr( T, 'plot','info_normest_norm_rho', 'v',verbose-1,'maxiteration',6 );
    fprintf( '|' ); tjsr( T, 'plot','polytope', 'maxiteration',5, 'v',verbose-1);
    %
    clf;
    T = tgallery( 'mejstrik_119', 20, 2 );
    fprintf( '|' ); tjsr( T, 'plot','tree', 'v',verbose-1, 'maxiteration',5, 'maxsmpdepth',8 );
    fprintf( '|' ); tjsr( T, 'plot','polytope', 'maxiteration',5, 'v',verbose-1, 'maxsmpdepth',8 );
    fprintf( '|' ); tjsr( T, 'plot','o', 'maxiteration',5, 'v',verbose-1, 'maxsmpdepth',8 );
    %
    % stupid options test
    fprintf( '|' ); 
    T = {[1 -1;1 0],[-1 1;0 -1]};
    tjsr( T, 'v',verbose-1, 'naturalselection',1, 'maxiteration',3, 'maxsmpdepth',1 );
    tjsr( T, 'v',verbose-1, 'maxnumcandidate',1, 'maxiteration',3 );
    %
    fprintf( '\n');
end



%% test tjsr_ 
%these files are not tested

%% test tliftproduct
EXPECT_EQ( unique( cell2mat(tliftproduct({2 3 5}, 2)) ), [1 2 3 4 5 6 9 10 15 25] );

%% test tliftsemidefinite
val = tliftsemidefinite({[1 1;0 1] [1 0; 1 0]}, 2);
EXPECT_EQ( val, {[1 4 4 2 4 1;0 1 2 1 3 1;0 0 1 0 2 1;0 0 0 1 2 1;0 0 0 0 1 1;0 0 0 0 0 1]; [ones(6,1) zeros(6,5)]} );

%% test trho
EXPECT_NO_THROW( 'trho([1 1;1 1]);' );
if( INIT('all') )
    MM = tgallery( 'rand_gauss', 10, 2, 'rho' );
    EXPECT_NEAR( trho(MM{1}), 1, 1e-10 );
    EXPECT_NO_THROW( 'trho(sparse(randi(2,2)));' );
    EXPECT_EQ( [-inf 1 2 1],trho({[], [1], [1 1;1 1], sparse([1 0;0 0])}) );
    EXPECT_NO_THROW( 'trho([]);' );
    EXPECT_NO_THROW( 'trho({1, [1 1;0 1]});' );
end

