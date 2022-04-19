% test subdivision

%#ok<*NASGU>
%#ok<*CTPCT>
%#ok<*ASGLU>

% preconditions


%% test version
EXPECT_FALSE( verLessThan( 'matlab','9.1') ); 

%% test toolboxes    
lic = license( 'test', 'symbolic_toolbox' );
if( ~lic );
    warning( 'testsubdivision:symbolic', 'ERROR: Symbolic Toolbox is not installed or licensed. The t-packages may not run.' ); end;

lic = license( 'test', 'signal_toolbox' );
if( ~lic );
    warning( 'testsubdivision:signal', 'ERROR: Signal Processing Toolbox is not installed or licensed. The t-packages may not run.' ); end;

%% test blf
EXPECT_NO_THROW( 'S_ = getS( ''a'',1/2*[1 2 1]'', ''M'',2 );' );
EXPECT_NO_THROW( 'blf( S_, ''iteration'',4, ''verbose'',0 );' );
if( INIT('all') )
    S_ = getS('2_butterfly');
    blf( S_,'iteration',3, 'wavelet', 'v',0 ); 
    blf( S_,'iteration',3, 'diff',1, 'verbose',0 ); 

    S_ = getS('1_121');
    blf( S_,'plot',0,'verbose',0);

    [c,PM,xyzv,oo_] = blf( {[],1}, S_, 'start',[1 -1]', 'removezero',1, 'start',[0 0 0 1 0 0 0]', 'iteration',5, 'verbose',0 );
    title( 'blf -- There should be a hat function.' );
    EXPECT_EQ( PM,32, oo_,{[],[1]}, size(xyzv,1),2 );
    try; 
        Ssym = getS('1_121'); 
        Ssym{1}.c = sym(Ssym{1}.c);
        blf( Ssym, 'iteration',4, 'verbose',0 );
    catch; 
        vprintf( 'ERROR: blf() may be broken.\n\n', 'cpr','err' ); end;

    warning( 'off', 'blf:dimension' );
    EXPECT_ERROR( 'blf(S_,''start'',[1 2],''iteration'',3,''wavelet'',''v'',0);' );
    warning( 'on', 'blf:dimension' );

    EXPECT_WARNING( 'blf( [1 1 1], S_, ''iteration'', 4, ''v'',0 );', 'blf:oolength' );
end


%% test checktile
M1_ = [1 2; -2 -2];  D1_ = [0 1; 0 -1]; 
M2_ = [-1 1; -1 -2]; D2_ = [0 0 0; -2 -1 0]; 
M12_ = M2_*M1_; 
D12 = setplus( M2_*D1_, D2_ );
[f1,C1,T1] = checktile(M1_,D1_,'verbose',-1,'legacy',0);
EXPECT_EQ( size(C1),[6 6]);
EXPECT_EQ( size(T1),[2 6]);
if( INIT('all') )
    f1=checktile(M1_,D1_,'verbose',-1);
    f2=checktile(M2_,D1_,'verbose',-1);
    f12=checktile(M12_,D12,'verbose',-1);
    EXPECT_EQ( [f1 f2 f12],[true, true, false]);
end


%% test compresscoordinates
cc = round( 3*compresscoordinates([10 20 ; 40 50],[2 1; -1 1]) );
EXPECT_EQ( cc, [0 1 -1 0;0 1 2 3;30 120 60 150]);

%% test constructdigit    
EXPECT_EQ( constructdigit(2),[0 1]);
if( INIT('all') )
    EXPECT_NO_THROW( 'constructdigit( [2 0; 1 2], ''ZZ'',[0 0; 1 0]'',''sym'');' );
    EXPECT_NO_THROW( 'constructdigit( [2 0; 1 2], ''ZZ'',[0 0; 1 0]'',''classify'',''sym'',''verbose'',0);' );
    EXPECT_EQ( size(constructdigit( [2 0; 1 2], 'random',4)),[2,4]);
    EXPECT_EQ( constructdigit( [2 0; 1 2] ),[0 0 1 1;0 1 1 2]);
    EXPECT_PRED( @isempty, constructdigit( 0 ) );
    EXPECT_NO_THROW( 'constructdigit( 1,''random'',2);' );
    EXPECT_NO_THROW( 'constructdigit( [10 9;11 9],''classify'',''ZZ'',[0 0;1 2;10 15]'',''sym'');' );
    EXPECT_NO_THROW( 'constructdigit( [10 9;11 9]);' );
    EXPECT_NO_THROW( 'constructdigit( sym([10 9;11 9]) );' );
    EXPECT_ERROR( 'constructdigit( [10 9;11 9],''classify'');' );
    EXPECT_ERROR( 'constructdigit( [10 9;11 9],''classify'',''ZZ'',[0 0;1 2;10 15]'',''random'',3);' );
end

%% test constructOmega
S_ = getS('1_1331');
Om_ = constructOmega(S_,'legacy');
EXPECT_EQ( sort(Om_),[0 1 2]);
if( INIT('all') )
    Om_ = constructOmega(S_,'Omega',[4 3 2 1 0],'stable');
    EXPECT_EQ( Om_,[4 3 2 1 0]);
    EXPECT_NO_THROW( 'constructOmega(getS(''2_butterfly''),''plot'');' );
    EXPECT_NO_THROW( 'constructOmega(S_);' );
    EXPECT_ERROR( 'constructOmega(S_,''Om'',[1 1.5]);' );
    EXPECT_ERROR( 'constructOmega(S_,''Om'',[1;2]);' );
end


%% test constructordering
oo_=constructordering([1 2],[3 ],[10 20],[30]);
EXPECT_EQ( oo_,{[1 2;10 20], [3;30]});
EXPECT_ERROR( 'constructordering(oo_,[1 2]);' );
EXPECT_PRED( @isordering, constructordering(2,2,'random',6) );


%constructPhi %untested

%constructU   %untested

%% test constructV/constructVt
Om_=[0 1 2 3];
V0_=constructV([1 1 1;0 1 1],'01');
EXPECT_PRED( @iscell, V0_ );
V0_=constructVt([1 1 1;0 1 1],'01');
EXPECT_PRED( @iscell, V0_ );
if( INIT('all') )
    EXPECT_EQ(constructV([],3),[]);
    EXPECT_EQ(constructV([],3),[]);
    EXPECT_EQ(constructVt([]),{[]});
    EXPECT_EQ(constructVt([]),{[]});

    EXPECT_EQ(constructV(sym(Om_)),constructV(Om_));
    EXPECT_EQ(constructV(sym(Om_)),constructV(Om_));
    EXPECT_EQ(constructVt(sym(Om_),0),constructVt(Om_,0));
    EXPECT_EQ(constructVt(sym(Om_),0),constructVt(Om_,0));
    EXPECT_PRED( @isempty, constructV(Om_,5) );
    [V,Om_]=constructVt([1 1 0 1 1],2,'01');
    EXPECT_PRED( @isempty, V );
    EXPECT_EQ( Om_,[0 0 0 0;0 1 3 4]);
    [V,Om_]=constructVt([1 1 0 1 1],2,'01');
    EXPECT_PRED( @isempty, V );
    EXPECT_EQ( Om_,[0 0 0 0;0 1 3 4]);
    [V,Om_]=constructVt([1 1 0 1 1],2,'01');
    EXPECT_PRED( @isempty, V );
    EXPECT_EQ( Om_,[0 0 0 0;0 1 3 4]);
    Om_=[0 0;2 2; 0 1;1 2;0 2]';
    Vt=constructVt(Om_,1); 
    V=constructV(Om_,1);
    EXPECT_EQ( rank(intersectspace(V,Vt)), 2 );
end

%% test daubechiesmask
EXPECT_PRED( @issym, daubechiesmask( 5 ) );
EXPECT_EQ( daubechiesmask( 0 ), 2 );

%differencescheme %untested

%% test dimVVt
dim_ = dimVVt([1 1 0 0;1 0 1 1;0 1 0 0],'01','verbose',0);
EXPECT_EQ(dim_,[5 3; 3 0]);
if( INIT('all') )
    Om_=[0 0;2 2; 0 1;1 2;0 2;1 0;2 0]';
    dim_=dimVVt(Om_,'V','verbose',0);
    dimt_=dimVVt(Om_,'Vt','verbose',0);
    EXPECT_EQ( dimt_,[6 3 ]);
    EXPECT_EQ( dim_,[6 4 1]);
end

%% test findperiod
[of_] = findperiod( [1 2 1 2 1 2; 4 5 5 5 5 5] );
if( INIT('all') )
    EXPECT_PRED( @isempty, findperiod(zeros(0,5)) );
    [oo_,pp_,nrep_]=findperiod([]);
    EXPECT_EQ( oo_,[], pp_,[], nrep_,0 );
    EXPECT_NO_THROW( 'findperiod([1]);' );
    EXPECT_NO_THROW( 'findperiod([1 2]);' );
    EXPECT_NO_THROW( 'findperiod(randi(2,1,500) );' );
    [oo_,pp_,nrep_]=findperiod([1 2 1 2 1 2; 4 5 5 5 5 5]);
    EXPECT_EQ( oo_,[1;4], pp_,[2 1; 5 5] );
    EXPECT_GE( nrep_, 1 );
    EXPECT_LE( nrep_, 3 );
    [of_]=findperiod([1 2 1 2 3 2 1 2 1 2],'fuzzy',1);
    EXPECT_PRED( @isempty, of_{1} );
    EXPECT_EQ( of_{2},[1 2] );
    [of_]=findperiod([1 2 1 2 1 2 1 2 2 2 ],'everywhere','verbose',0);
            EXPECT_PRED( @isempty, of_{1} );
    EXPECT_EQ( of_{2},[1 2] );
end


%% test getS
EXPECT_NO_THROW( 'getS(''a'',1/2*[1 2 1]'',''M'',2);' );
if( INIT('all') )
    EXPECT_WARNING( 'getS( ''a'',1/2*[1 2 1],''M'',2)', 'getS:transpose' );
    EXPECT_WARNING( 'getS( ''a'',[1 2 1]'',''M'',2);', 'getS:sumrule0' );
    EXPECT_WARNING( 'getS( ''a'',1/2*[1 2 1]'',''M'',2,''D'',[0 2],''bigcheck'')', 'getS:digit', 'getS:bigcheck' );
    EXPECT_WARNING( 'getS( ''a'',1/2*[1 1]'',''M'',1,''D'',[0 2],''bigcheck'');', 'getS:notexpanding', 'getS:bigcheck' );
    EXPECT_WARNING( 'getS( ''a'',1/2*[1 1]'',''D'',[0 2]);', 'getS:M' );

    warning( 'off', 'normalizeS:support' );
    EXPECT_NO_THROW( 'getS(''a'',[1 2 1],''M'',2,''name'',''tommsch'',''nocheck'');' );
    S1_ = getS('a',1/2*[1 2 1]','M',2);
    EXPECT_EQ( S1_,getS(S1_) );
    EXPECT_EQ( S1_,getS({1/2*[1 2 1]',2}) );
    EXPECT_NO_THROW( 'getS(''1_121'');' );
    EXPECT_NO_THROW( 'getS(S1_,''supp'');' );
    EXPECT_NO_THROW( 'getS([S1_; S1_],''supp'');' );
    EXPECT_NO_THROW( 'getS(S1_,''OmegaRR'');' );
    EXPECT_NO_THROW( 'getS(S1_,''characteristic'');' );
    EXPECT_NO_THROW( 'getS(''1_rand'');', 'normalizeS:support' );
    EXPECT_NO_THROW( 'getS(''2_rand'');', 'normalizeS:support'  );
    EXPECT_NO_THROW( 'getS(''1_all'',''nocheck'');', 'normalizeS:support' ); %1_all also includes non-stationary schemes %XX
    S2_ = getS('2_rand');
    warning( 'on', 'normalizeS:support' );

    EXPECT_WARNING( 'getS([S1_;S2_]);', 'getS:dimM' );
    S2_ = getS('1_1133');
    S_ = getS([S1_;S2_]);

    EXPECT_WARNING( 'getS(100,''nocheck'')', 'getS:nodata' );

    val_ = getS('1_143'); 
    EXPECT_DOUBLE_EQ( val_{1}.c, 1/4*[1;4;3] );
    EXPECT_NO_THROW( 'getS(''3_rand'');', {'normalizeS:support'} );
    EXPECT_NO_THROW( 'getS(val_);' );
    EXPECT_NO_THROW( 'getS(''a'',1/2*[1 2 1]'',''M'',2); ' );
    EXPECT_EQ( val_{3},[0 1]);
    EXPECT_NO_THROW( 'getS(S1_);' );
end


%% test isordering
S1_ = getS('1_121');
EXPECT_EQ( isordering(findperiod([1 2 1 2 1 2; 4 5 5 5 5 5])), true);
EXPECT_EQ( isordering([2 3]),false);
EXPECT_EQ( isordering([]),false);
EXPECT_EQ( isordering({[1 3 2]}),false);



%% test isS
S1_ = getS('1_121');
EXPECT_EQ( isS({[],[1 2],[]}),false);
EXPECT_PRED( @isS, S1_ );
EXPECT_NPRED( @isS, [2 3] );
EXPECT_PRED( @isS, cell(0,4) );


%% test isT
T_ = transitionmatrix( getS('1_1133') );
EXPECT_PRED( @isT, T_ );


%% test multiplyS
S1_ = getS('1_143');
S2_ = getS('1_1133');
EXPECT_WARNING( 'multiplyS();', 'getS:notexpanding' );
EXPECT_NO_THROW( 'multiplyS(S1_); ' );
EXPECT_NO_THROW( 'multiplyS(S1_,S2_); ' );
EXPECT_NO_THROW( 'multiplyS(S1_,S2_,S2_);' );
if( INIT('all') )
    EXPECT_ERROR( 'multiplyS(2,3)' );
    S_ = getS('2_not_jointly_expanding');
    EXPECT_WARNING( 'multiplyS(S_(2,:),S_(1,:),S_(2,:))' ,'getS:notexpanding', 'multiplyS:notexpanding' );
end


%% test normalizeS
S_ = getS({[1 1 1]',2},'nocheck');
EXPECT_NO_THROW( 'normalizeS(S_,''scale'');' );
if( INIT('all') )
    EXPECT_WARNING( 'normalizeS({[1 0 1],2});', 'normalizeS:support' );
    EXPECT_NO_THROW( 'normalizeS(S_,''scale'');' );
    EXPECT_NO_THROW( 'normalizeS(S_,''gauss'');' );
    Sn_ = normalizeS(S_,'equal');
    EXPECT_EQ( Sn_{1,1}.c,1/2*[1 2 1]');
    EXPECT_WARNING( 'normalizeS(getS({[1 2 -1]'',2},''nocheck''));', 'normalizeS:failure' );
end


%% test num2ordering
S_=getS({[1 1]',2},'nocheck');
Snum_=getS({[],[1 2; -1 2]});
pt_=[0.1875;0.375];
EXPECT_NO_THROW( 'ordering2num(num2ordering(Snum_,pt_,''l'',50),Snum_);' );
if( INIT('all') )
    [do_, doint_, num_, proofnum_] = num2ordering([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],getS('M',2),.4,'proof','v',0);
    EXPECT_DOUBLE_EQ( proofnum_, .4 );

    num_=ordering2num(num2ordering(Snum_,pt_,'l',50,'verbose',0,'plot',1),Snum_);
    EXPECT_DOUBLE_EQ( num_, pt_ );
    EXPECT_WARNING( 'num2ordering(Snum_,[1 2],''check'',1)', 'num2ordering:pointoutside' );

    [do_,doint_,num_,proofnum_]=num2ordering(Snum_,pt_,'check',1,'proof','verbose',0);
    EXPECT_EQ(do_,{[1 1;2 3],[1 1;1 2]});
    EXPECT_EQ(proofnum_,num_);
    Sn_={[],2,[0 1]; [],2,[0 3]};

    EXPECT_NO_THROW( 'ordering2num(num2ordering({[],[1 2]},Sn_,5/3,''plot'',1,''verbose'',0,''l'',10,''delta'',.7,''round'',0),Sn_);' );

    EXPECT_WARNING( 'num2ordering([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],Snum_,[.4 .2],''v'',0);', 'num2ordering:period' );
    EXPECT_ERROR( 'num2ordering([S_;S_],[.4 .2],''v'',0);' );
    EXPECT_WARNING( 'num2ordering(getS(''M'',10),0.1133333333333333333333333333,''l'',4);', 'num2ordering:period' );
end

%% test ordering2num
num_=ordering2num({[2 1 ],[]},getS({[],2}) );
EXPECT_EQ(num_,1/2);
if( INIT('all') )
    EXPECT_DOUBLE_EQ( ordering2num({[1 2],[1 2]},10), 0.010101010101010101 );
    EXPECT_NO_THROW( 'ordering2num({[1 2],[1 2]},[1 2;2 2]);' );
    EXPECT_NO_THROW( 'ordering2num({[2 1],[1 3]},getS({[],[1 2; -1 1]}) );' );
    num_ = ordering2num({[],[1 2]},getS({[],2}) );
    EXPECT_EQ(num_,1/3);
    num_=ordering2num({[],[]},getS({[],2}) );
    EXPECT_EQ(num_,0);
    S_=getS({[],2});
    EXPECT_ERROR( 'ordering2num({[2 1 ]},S_);' );
    EXPECT_ERROR( 'ordering2num({[1],[2 1 ]},[S_;S_]);' );
end

%% test ordering2vector
vec=ordering2vector({ [10 20 30] , [1 2] },10);
EXPECT_EQ( vec,[10 20 30 1 2 1 2 1 2 1]);
if( INIT('all') )
    EXPECT_NO_THROW( 'ordering2vector({[],[1]},5);' );
    vec=ordering2vector({[1 2 3 4 5],[]},10);
    EXPECT_EQ(vec,[1 2 3 4 5]);
    EXPECT_NO_THROW( 'ordering2vector({[],[]});' );
end

%% test peter
EXPECT_NO_THROW( 'peter(rand(2,30),4);' );
EXPECT_NO_THROW( 'peter(rand(2,30),.4);' );


%% test restrictmatrix
S_ = getS( '1_121' );
[T_,Om_,V0_] = transitionmatrix( S_, 'V',0, 'Omega',[0 1 2 3] );
EXPECT_NO_THROW( 'restrictmatrix( T_{1}, V0_ );' );
if( INIT('all') )
    S_ = getS('1_121');
    [T_,Om_,V0_] = transitionmatrix( S_, 'V',0, 'Omega',[0 1 2 3] );
    EXPECT_NO_THROW( 'restrictmatrix( T_{1}, V0_ );' );
    EXPECT_WARNING( 'restrictmatrix(T_,[1 -3 3 1]'');', 'restrictmatrix:notinvariant', 'restrictmatrix:notest' );

    [TA,TT,TR,NULL,BASIS] = restrictmatrix(T_,[1 -1 0 0; 1 0 -1 0; 1 -1 0 0]','smallsize',1 );
    syms xxx_; 
    EXPECT_ERROR( 'restrictmatrix(sym([1 2;-1 -2]),[xxx_ -1]'');' );
end

%testsubdivision %this file

%% test tile
%figure; 
S_=getS('2_butterfly');
[Q,oo_] = tile(S_,'verbose',0,'peter',1000,'round',[1e-1 1e-5]); 
EXPECT_EQ( size(Q,1), 2);
EXPECT_FALSE( any(diff(oo_)) );
if( INIT('all') )
    EXPECT_NO_THROW( 'Q=tile(S_,''verbose'',0,''peter'',[.4 .8]); ' );
    EXPECT_NO_THROW( '[Q,oo_]=tile(S_,''verbose'',0,''plot'',{''MarkerSize'',10});' );
    EXPECT_NO_THROW( 'tile(S_,''verbose'',0,''digit'');' );
    EXPECT_NO_THROW( 'tile(S_,''verbose'',0,''diffdigit'',''iteration'',3); ' );
    EXPECT_NO_THROW( 'tile(S_,''verbose'',0,''supp''); ' );
    EXPECT_NO_THROW( 'tile(S_,''verbose'',0,''OmegaRR'',''iteration'',4); ' );
    EXPECT_WARNING( 'tile({[1 2 2 1],[]},[getS(''2_rand'');S_],''verbose'',0,''iteration'',5);','tile:oolength' );
    tile([getS('2_rand');getS('2_rand')],'supertile',6,'round',[1e-2 1e-4],'verbose',0); 
    tile([getS('2_rand'); getS('2_rand')],'supertile',6,'iteration',8,'verbose',0); title('tile - supertile (something Rohrschach-like)');
    tile(getS('M',3,'D',[0 1 5]),'interval','start',{[0 1]},'round',0,'v',-1); title('tile - simple interval arithmetic');
end

%% test tilearea
S_=getS('2_butterfly');
EXPECT_EQ( tilearea(tile(S_,'plot',0,'verbose',0),'verbose',0), 1);
if( INIT('all') )
    S_=getS('2_butterfly');
    tilearea(S_,'verbose',0);
    tilearea(getS('M',3,'D',[0 1 5]),'v',0,'visual');        
    tilearea(tile(S_,'plot',0,'v',0),'v',0,'visual');
end

%% test transitionmatrix
S1_=getS('1_121');
S2_=getS('1_1133');
S_=[S1_;S2_];
EXPECT_NO_THROW( '[T_,Om_]=transitionmatrix(S_);' );
if( INIT('all') )
    transitionmatrix(S_,'colsum',0);
    transitionmatrix(S_,'colsum',1);
    EXPECT_WARNING( 'Om_ = transitionmatrix(S_,''colsum'',2); ', 'transitionmatrix:colsum' );
    EXPECT_EQ( Om_,[0 1 2]);
    EXPECT_ERROR( 'transitionmatrix(S_,''Omega'',[0;0]);' );

    EXPECT_WARNING( 'transitionmatrix(S_,''Omega'',[50],''colsum'',1);', 'transitionmatrix:colsum' );

    T_=transitionmatrix( S_, 'noflat' );
    EXPECT_PRED( @iscell, T_{1} );

    [T_,Om_,V0_]=transitionmatrix( S_, 'V',0 );
    EXPECT_EQ( null(V0_','r'), [1;1;1] );

    [T_,Om_,V1]=transitionmatrix( S_, 'V',3 );
    EXPECT_EQ(Om_,[0 1 2 3 4]);
    EXPECT_NO_THROW( 'transitionmatrix(S_,''noflat'',''colsum'',1);' );
    EXPECT_NO_THROW( 'transitionmatrix(S_,''onlyindex'');' );
    EXPECT_NO_THROW( 'transitionmatrix(S_);' );
end

%% test vector2ordering
[o1, o2] = vector2ordering( [1 2 1 2 1 2; 4 5 5 5 5 5] );
EXPECT_EQ( o1,[1;4], o2,[2 1; 5 5] );
