% test tmisc

%#ok<*NASGU>
%#ok<*CTPCT>
%#ok<*ASGLU>

syms zz;
syms xx; 
assume( xx, 'real' );

% preconditions

%% test flatten
EXPECT_EQ( flatten({{2 3 {4; 5} 6} {{{7}}}}), {2 3 4 5 6 7} );
flatten({2});
EXPECT_PRED( @iscell, flatten([2 3]) );

%% test grCenter
G_ = [1 2;2 5; 5 6; 2 3; 3 4; 4 7; 4 8; 8 9];
EXPECT_EQ( grCenter(G_), 3 );
EXPECT_EQ( grCenter(G_,'edge'), [3;2] );
G_ = [1 2; 2 3; 4 5; 2 4];    
EXPECT_EQ( sort(grCenter(G_)), [2 4] );


%% test grVerCover
G_ = [1 2; 1 3; 1 4; 2 5; 1 6];
G_ = digraph( G_(:,1), G_(:,2) );
val_ = grVerCover(G_);
EXPECT_FALSE ( numel(val_)~=3 || ~searchincellarray([1 2],val_,1) || ~searchincellarray([1 5],val_,1) || ~searchincellarray([2 3 4 6],val_,1) );
G_ = [1 1; 1 3; 1 2; 2 1];
G_ = digraph( G_(:,1), G_(:,2) );
val_ = grVerCover( G_ );
EXPECT_EQ( val_, {1} );


%% test identifymatrix
EXPECT_NO_THROW( 'identifymatrix( {[1 2 3; 0 2 1; 0 0 1],[2]}, 2 );' ); %#ok<NBRAK>
EXPECT_EQ( identifymatrix(1i,2,'nan'), 0 );
if( INIT('all') && exist('tgallery','file') )
    EXPECT_FALSE( identifymatrix(xx,'pm1') );
    EXPECT_FALSE( identifymatrix(zz,'bool') );
    EXPECT_TRUE( identifymatrix(1,'pm1') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_stochastic',4,2),'columnstochastic') );


    EXPECT_TRUE( identifymatrix(tgallery('rand_stochastic',4,2),'columnstochastic') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_stochastic',4,1,'nocell')','rowstochastic') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_doublestochastic',4,2),'doublestochastic') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_pm1',4,2),'pm1') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_bool',4,2),'bool') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_hess',4,2),'hessu') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_hess',4,1,'nocell')','hessl') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_unitary',4,1),1,'unitary') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_zero',4,1),1,'zero') );
    EXPECT_TRUE( identifymatrix(tgallery('rand_unitary',4,1),1,'inv') );
    EXPECT_TRUE( identifymatrix(tgallery('orthog_2',4,1),1,'unitary') );
    EXPECT_TRUE( identifymatrix(sym([1 2;2 1]),'sym') );

    EXPECT_WARNING( 'identifymatrix([1 2;3 4],''wrongname'')',  'identifymatrix:wrongname' );
end
%XX more tests for identifymatrix


%% test intersectspace
[in_,dim_] = intersectspace( [1 1 1; 1 2 1; 1 2 2]', [1 1 1; 1 2 1]', [1 1 1]' );
EXPECT_EQ( in_, [1;1;1] );
EXPECT_EQ( dim_, 1 );
[in_,dim_] = intersectspace( [1 -1 0 0; 1 0 -2 1; 0 1 -2 -1]', [1 1 0 0; 1 0 -1 0]' );
EXPECT_EQ( rank(in_), 1 );
EXPECT_EQ( dim_, 1 );
[in_,dim_] = intersectspace([1 0 0; 0 1 0;1 1 0]');
EXPECT_EQ( dim_, 2 );
[in_,dim_]=intersectspace([1 0 0]',[0 1 0;1 -1 0]',[0 1 1;0 1 -1]',[1 1 1;1 1 -1; 1 -1 1]');
EXPECT_EQ( dim_, 0 );
[in_,dim_]=intersectspace(sym([1 0 0]'),[1 1 0; 1 -1 0]');
EXPECT_EQ( dim_, 1 );
if( INIT('all') )
    EXPECT_ERROR( 'intersectspace()' );
    EXPECT_ERROR( 'intersectspace([2 3 3;1 1 1]'',[1 2 ]'');' );
end;


%% test issquare
EXPECT_PRED( @issquare, rand(3,3,3) );
EXPECT_NPRED( @issquare, rand(3,4)  );


%% test issym
EXPECT_PRED( @issym, sym('23.2') );
EXPECT_NPRED( @issym, 3.2 );


%% test iswholenumber
EXPECT_EQ( iswholenumber([2 2.1 inf 0]), [true false false true] );
EXPECT_EQ( iswholenumber(@sum), false );
EXPECT_EQ( iswholenumber(sym(2)), true );
EXPECT_EQ( iswholenumber(sym(2.5)), false );
EXPECT_EQ( iswholenumber(xx), false );


%% test lexicographic
EXPECT_EQ( lexicographic( [0 1 0; 0 1 2; 3 1 1] ,2), [1 0 0;1 2 0;1 1 3] );
EXPECT_NO_THROW( 'lexicographic( [0 1 0; 0 1 2; 3 1 1]);' );
EXPECT_NO_THROW( 'lexicographic( [0 1 0; 0 1 2; 3 1 1],''inf'');' );
EXPECT_ERROR( 'lexicographic( [0 1 0; 0 1 2; 3 1 1] ,@(x) x^2)' );


%% test liminf
EXPECT_EQ( liminf([10 1 9 2 8 3 7 4 6]), [1 1 2 2 3 3 4 4 6] );

%% test limsup
EXPECT_EQ( limsup([10 1 9 2 8 3 7 4 6]), [10 9 9 8 8 7 7 6 6] );


%% test makepositive
EXPECT_EQ( makepositive([ 0 -1 2]), [0 1 -2] );
EXPECT_EQ( makepositive(sym([ 0 -1 2])), sym([0 1 -2]) );
EXPECT_DOUBLE_EQ( makepositive(1i), 1 );


%% test mergeinterval
[lower_, upper_] = mergeinterval( [0 1 2 3 4], [1.5 1.6 3.5 3 5] );
EXPECT_EQ( lower_, [0 2 4] );
EXPECT_EQ( upper_, [1.6 3.5 5] );
EXPECT_NO_THROW( 'mergeinterval(sym([0 1.5 1 1.6 2 3.5 3 3 4 5]));' );
EXPECT_NO_THROW( 'mergeinterval({[0 1.5];[1 1.6];[2 3.5];[3 3];[4 5]});' );
EXPECT_NO_THROW( 'mergeinterval([0 1.5; 1 1.6; 2 3.5; 3 3; 4 5]);' );
EXPECT_NO_THROW( 'mergeinterval([0 1.5; 1 1.6; 2 3.5; 3 3; 4 5],''output'',''C'');' );
EXPECT_NO_THROW( 'mergeinterval({[]});' );
EXPECT_NO_THROW( 'mergeinterval([]);' );
EXPECT_NO_THROW( 'mergeinterval([]);' );
[~,~] = mergeinterval(); 
EXPECT_EQ( mergeinterval([0 2; 1 -1],'output','M'), [0 2] );
EXPECT_EQ( mergeinterval([0 2; -1 1],'output','M','switch'), [-1 2] );
EXPECT_EQ( mergeinterval([-inf inf],'output','M','switch'), [-inf inf] );
EXPECT_EQ( mergeinterval([],'output','M','switch'), zeros(0,2) );
if( INIT('all') )
    EXPECT_WARNING( 'mergeinterval([0 1 2 3 4],[1.5 1.6 3.5 3 5]);', 'mergeinterval:nargout' );
    EXPECT_ERROR( 'mergeinterval([0 1 2 3],[1.5 1.6 3.5 3 5]);' );
    EXPECT_ERROR( 'mergeinterval([0 1 2 3],[1.5 1.6 3.5 3 5]);' );
    EXPECT_ERROR( 'mergeinterval([0 1 2 3 3]);' );
    EXPECT_ERROR( 'mergeinterval([0 1 2; 3 3 3]);' );
    EXPECT_ERROR( 'mergeinterval([0 1; 3 3],''o'',''Wrongstring'');' );
    EXPECT_ERROR( 'mergeinterval([1 2 3;4 5 6;7 8 9]);' );
end

%% test mixvector
EXPECT_EQ( sortrows(mixvector(0:1, 2)')', sortrows([0 0 1 1;0 1 0 1]')' );
EXPECT_EQ( mixvector([-1 1], 3,0), 8 );
EXPECT_EQ( mixvector([-1 1], 3,1), [-1;-1;-1] );
[~,val_] = mixvector(0:1, 2);
EXPECT_EQ( val_, 4 );

EXPECT_EQ( mixvector([10 20 30 40],1), [10 20 30 40] );
EXPECT_EQ( mixvector([10 20 30 40],1,1), 10 );
EXPECT_EQ( mixvector([10 20 30 40],1,0), 4 );
if( INIT('all') )
    EXPECT_ERROR( 'mixvector([10 20 30 40],1,5);' );
    EXPECT_EQ( mixvector([1 2 3],2,5), [2;2] );
end


%% test nestedcellfun
EXPECT_EQ( nestedcellfun(@(x) sum(x),{[1 2],[2 3 4]}), [3 9] );
EXPECT_EQ( nestedcellfun(@(x) ndimsm(x),{{[1 2 1]},[1 ;2 ]},'UniformOutput',false), {{2},1} );


%% test nondiag
val_ = nondiag({[1 2] [2 3] ;[1] [40] },1); %#ok<NBRAK>
EXPECT_EQ( val_, {[1 2] []; 1 40} );
nondiag([1 2 3; 3 2 3; 3 4 5]);
EXPECT_NO_THROW( 'nondiag([1 2 3; 3 2 3; 3 4 5],-inf);' );
EXPECT_NO_THROW( 'nondiag([1 2 3; 3 2 3; 3 4 5],inf);' );


%% test normalizematrix
MM_ = [1 2 3; 0 -1 2];
val_ = normalizematrix(MM_,'colsum');
EXPECT_DOUBLE_EQ( sum(val_,1), 1 );
if( INIT('all') )
    EXPECT_NO_THROW( 'normalizematrix( [], ''colsum'' );' );
    val_ = normalizematrix( MM_, 'rowsum' );
    EXPECT_DOUBLE_EQ( sum(val_,2), 1 );
    EXPECT_NO_THROW( 'normalizematrix( MM_, ''dirsum'',3 );' );
    EXPECT_DOUBLE_EQ( normalizematrix(MM_,'dirmax',1), [1 1 1;0 -.5 2/3] );

    val_ = normalizematrix( MM_, 'dirnorm',[2 2] );
    EXPECT_DOUBLE_EQ( sum(val_.^2,2), [1;1] );

    MM_ = [-1 2; 0 3];
    val_ = normalizematrix(MM_,'rho');
    EXPECT_DOUBLE_EQ( rho(val_), 1 );

    val_=normalizematrix(MM_,'norm',inf);
    EXPECT_DOUBLE_EQ( norm(val_,inf), 1 );
    EXPECT_EQ( normalizematrix(MM_,'binary'), [1 1;0 1] );
    EXPECT_EQ( normalizematrix(MM_,'positive',2),[1 -2;0 3] );
    EXPECT_EQ( normalizematrix({[1 2],[-4 2]},'dotprod',{[1 1],[0 1]}), {[1/3 2/3],[-2 1]} );
    EXPECT_EQ( normalizematrix([1 2],'dotprod',[1 1]), [1/3 2/3] );
    EXPECT_ERROR( 'normalizematrix([10 20 30 40],''wrongstring'');' );
end


%% test num2color
EXPECT_NO_THROW( 'num2color(0); ');
EXPECT_EQ( num2color(4), num2color(4) );
EXPECT_NO_THROW( 'num2color(200,100);' );


%% test removezero
AA_ = [0 0 0 0;0 1 0 0;0 0 0 0;0 2 0 3];
EXPECT_EQ( removezero(AA_,'compact',1), [0 1 0 3;0 2 0 0] );
EXPECT_EQ( removezero(AA_,'compact',2), [0 0;1 0;0 0;2 3] );
EXPECT_EQ( removezero(AA_,'compact',[1 2]), [1 3; 2 0] );
EXPECT_EQ( removezero(AA_,[1 2]), [1 0; 2  3] );
EXPECT_EQ( removezero(AA_,[]), AA_ );
EXPECT_EQ( removezero(AA_,'border'), [1 0 0;0 0 0;2 0 3] );
EXPECT_EQ( removezero(AA_,'outside'), [0 0 0 0; 0 1 0 0;0 0 0 0;0 2 0 3] );
EXPECT_EQ( removezero(AA_,'inside'), [1 0 0;0 0 0;2 0 3] );
EXPECT_EQ( removezero(AA_,'all'), [1;2;3] );
EXPECT_EQ( removezero(AA_,'left'), [0 0 0; 1 0 0;0 0 0;2 0 3] );
EXPECT_EQ( removezero(AA_,'right'), AA_ );
EXPECT_EQ( removezero(AA_,'top'), [0 1 0 0;0 0 0 0;0 2 0 3] );
EXPECT_EQ( removezero(AA_,'bottom'), AA_ );
AA_ = [0 1 0; 0 2 0];
[val_,idx_]=removezero(AA_,[2],'keepdim'); %#ok<NBRAK>
EXPECT_EQ( val_, [1 0;2 0] );
EXPECT_EQ( idx_, [1;2] );
EXPECT_ERROR( 'removezero(AA,''wrongstring'');' );
AA_ = [0 0 1 0 1];
[~,idx_] = removezero( AA_, 'border' );
EXPECT_EQ( idx_, [1;3] );


%% test repcell
EXPECT_EQ( repcell(10 ,2), {10 10;10 10} );


%% test savetocellarray
val_ = savetocellarray( [1 2 3], logical([1 1 0 0 1 ]), {[-1 -2], [-3 -4 -5]} );
EXPECT_EQ( val_, {[1 2], [-3 -4 3]} );


%% test searchincellarray
[found_, ioo_] = searchincellarray( [1 2 3 4], {}, 0 );
EXPECT_EQ( found_, 0 );
EXPECT_EQ( ioo_, 0 );
[found_, ioo_] = searchincellarray( [1 2 3 4], {[1 2; 3 4],[1 2 3 4],[2 3 1]',[5]}, 0 ); %#ok<NBRAK>
EXPECT_EQ( found_, true );
EXPECT_EQ( ioo_, 2 );


%% test setplus
EXPECT_NO_THROW( 'setplus( [1 0; 1 1], [10 10; 20 20; 30 30], ''rows'' );' );
EXPECT_NO_THROW( 'setplus(1);' );
val_ = setplus( [1 0; 1 1]', [10 10; 20 20; 30 30]' );
EXPECT_EQ( val_, [11 11 21 21 31 31;10 11 20 21 30 31] );
val_ = setplus( [1 0], [0 1], 'stable' );
EXPECT_EQ( val_, [1 0 2] );
val_ = setplus( [1 0], [0 1], [0], [1], 'nounique' ); %#ok<NBRAK>
EXPECT_EQ( sort(val_), [1 2 2 3] );


%% test smallestchoice
EXPECT_NO_THROW( 'smallestchoice({[1 3;4 2], [2;3], [1;0], [2 3;6 1]});' );
EXPECT_NO_THROW( 'smallestchoice({});' );
EXPECT_NO_THROW( 'smallestchoice({randi(30,2,3), randi(30,2,3), randi(30,2,3)});' );


%% test subsco
val_ = subsco( [1 2 3; 4 5 6], [1 2; 2 3]' );
EXPECT_EQ( val_, [2 6] );
val_ = subsco( [1 2 3; 4 5 6], [1 2; 2 3]', [-1 -1] );
EXPECT_EQ( val_, [1 -1 3;4 5 -1] );
EXPECT_ERROR( 'subsco( [1 2 3; 4 5 6], [1 2; 4 3]'' );' );
val_ = subsco([1 2 3; 4 5 6],[1 2; 4 3]','save');
EXPECT_ERROR( 'val = subsco( [1 2 3; 4 5 6], [1 2; 1 2]'', [1 2], ''save'' );' );


%% test tbuildproduct
MAT_ = {[0.8 0.2; 0.3 0.7],[0.4 0.6;0.5 0.5]};
EXPECT_NO_THROW( 'tbuildproduct(MAT_, {[],[1 2]});' );
if( INIT('all') )
    SOL_ = 1/21*[11 10;11 10];
    EXPECT_DOUBLE_EQ( tbuildproduct(MAT_,{[],[1 2]}), SOL_ );
    EXPECT_DOUBLE_EQ( tbuildproduct(MAT_,{[],[1 2]},'sym'), SOL_ );
    EXPECT_DOUBLE_EQ( tbuildproduct(MAT_,{[],[2 1]},'reverse'), SOL_ );
    EXPECT_NO_THROW( 'tbuildproduct( MAT_, {[],[1 2]}, ''l'',1 );' );

    SOL_ = 1/20*[10 10;11 9];
    EXPECT_DOUBLE_EQ( tbuildproduct(MAT_, [1 2]), SOL_ );
    EXPECT_DOUBLE_EQ( tbuildproduct(MAT_, {[1 2],[]}), SOL_ );
    EXPECT_ERROR( 'tbuildproduct({{[10],[20]}},[1 2]);' );
    EXPECT_EQ( tbuildproduct({{[10],[20]}},[1 1;1 2]), 200 );
end


%% test tif
EXPECT_EQ( [tif(true, 1, 2) tif(false, 1, 2)], [1 2] );


%%test unflatten/flatten
C_ = {[2 3 4] [5 ; 5 ; 5 ] {3 7; 4 6}; 40 50 {60  80}};
EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_={1};              EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_=[1 2 3;4 5 6];    EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_={1 {} };          EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_={};               EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_={{}};             EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_=[];               EXPECT_EQ( C_, unflatten(flatten(C_),C_) );
C_={1 2 3};          EXPECT_EQ( unflatten( 10, C_, 'assign' ), {10 10 10} );
C_={};               EXPECT_EQ( unflatten( 10, C_, 'assign' ), {} );
C_={{}};             EXPECT_EQ( unflatten( 10, C_, 'assign' ), {{}} );


%% test uniquecell
val_ = uniquecell( {1 2 1 [1 2] [1 3] [1 2]} )';
EXPECT_EQ( numel(val_), 4 );
[AUC_, IDXC_ ,IDXCC_] = uniquecell({1 2 3 1 3 4 5});
[AUA_, IDXA_ ,IDXAA_] = unique(    [1 2 3 1 3 4 5]);
EXPECT_EQ( cell2mat(AUC_), AUA_ );
EXPECT_EQ( IDXC_, IDXA_ );
EXPECT_EQ( IDXCC_, IDXAA_ );


%% test vdisp
val_ = vdisp( {[1 2 3]} );
EXPECT_TRUE( ischar(val_) && strfind(val_,'1') && strfind(val_,'2') && strfind(val_,'3') ); %#ok<STRIFCND>
if( INIT('all') )
    ST_ = struct; 
    ST_.There = 'should'; 
    ST_.be = 'a'; 
    ST_.struct = @displayed; 
    val_ = vdisp(ST_);
    EXPECT_TRUE( contains(val_,'struct') && contains(val_,'@displayed') && contains(val_,'There') );
    val_ = vdisp({[123 234 345]},'OLD');
    EXPECT_TRUE( ischar(val_) && contains(val_,'OLD') && contains(val_,'123') && contains(val_,'234') && contains(val_,'345') );
end

%% test vprintf
val_ = vprintf('12345', 'cpr', [.8 .4 .24], 'str','', 'npr' );
EXPECT_TRUE( contains(val_,'12345') );
if( INIT('all') )
    vprintf('ERROR: vprintf() is broken.\n\n', 'imp',[3 2], 'cpr','err' );
    vprintf('ERROR: vprintf() is broken.\n\n', 'cpr','err', 'noprint' );
    val_ = vprintf('TEST','cpr','err','sze',[10 100], 'noprint' );
    EXPECT_FALSE( contains(val_, 'TEST') );
    val_ = vprintf('%v',[1 2 3;4 5 6], 'npr' );
    newl_ = sprintf('\n'); %#ok<SPRINTFN>
    EXPECT_TRUE( contains(val_,'1') && contains(val_,'6') && contains(val_,newl_) );
    stru_.num=6;
    val_ = vprintf( 'Numbers from 1 to 6 : \n   %i %f %v %r %v\n', 1, 2, {{3}}, [4 5], stru_, 'npr',1 );
    EXPECT_TRUE( contains(val_,'1') && contains(val_,'2.0') && contains(val_,'3') && contains(val_,'4') && contains(val_,'5') && contains(val_,'num') && contains(val_,'6') );
    val_ = vprintf('Some text again.', 'str','OLD', 'noprint' );
    EXPECT_EQ( val_, 'OLDSome text again.' );
    val_ = vprintf( '%i%%',100, 'noprint' );
    EXPECT_EQ( val_, '100%' );
    EXPECT_ERROR( 'vprintf(''ERROR: vprintf() is broken. Wrong number of arguments %i\n\n'',''cpr'',''err'');' );
end

%% test R2017b
try;
    val = runtests('testtmiscR2017b');
    if( ~all([val.Passed]) );
        assert( false ); end;
catch
    assert( false );
end;

