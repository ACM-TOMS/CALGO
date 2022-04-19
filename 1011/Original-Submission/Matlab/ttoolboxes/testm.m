% test sequence

%#ok<*CTPCT>
%#ok<*ASGLU>
%#ok<*NASGU>

if( exist('symbol2mask','file')~=0 ); 
    sequenceinstalled = 1; 
else; 
    fprintf( 'Warning: sequence-package not installed.\n' );        
    sequenceinstalled = 0;  end;

if( exist('vprintf','file')~=0 ); 
    tmiscinstalled = 1; 
else; 
    fprintf( 'Warning: tmisc-package not installed.\n' );        
    tmiscinstalled = 0; end;

if( opengl('info') );
    ogl = 1;
else;
    ogl = 0; end;



% preconditions    

%% test allm
AA = [1 1 2;0 1 1; 1 1 1];
EXPECT_EQ( allm(), 1);
EXPECT_EQ( allm(AA), 0 );
EXPECT_EQ( allm(AA,1), [0 1 1] );
EXPECT_EQ( allm(AA,2), [1;0;1] );
EXPECT_EQ( allm(AA,[1 2]), 0 );
EXPECT_EQ( allm([1 0;1 2],[]), [1 0;1 1] );

%% test anym
BB = [0 1 2;0 1 1; 0 1 1];
EXPECT_EQ( anym(), 0 );
EXPECT_EQ( anym(BB), 1 );
EXPECT_EQ( anym(BB,1), [0 1 1] );
EXPECT_EQ( anym(BB,2), [1;1;1] );
EXPECT_EQ( anym(BB,[1 2]), 1 );
EXPECT_EQ( anym([1 0;1 2],[]), [1 0;1 1] );


%% test sph2cartm/cart2sphm/sph2cartm2/cart2sphm2
EXPECT_DOUBLE_EQ( cart2sphm([1 0]'), [0 1]' );
EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 1 0]'), [pi/2 0]' );
val = randn(100,300); 
val = normalizematrix( val, 'colnorm', 2 ); 
EXPECT_NEAR( sph2cartm2(cart2sphm2(val)), val, 5e-12 );
if( INIT('all') )
    %S^1
    %cart2sphm, type=0
    EXPECT_DOUBLE_EQ( cart2sphm([1 0]'),[0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 1]'),[pi/2 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([-1 0]'),[pi 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 -1]'),[-pi/2 1]');
    %cart2sphm2, type=0
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[1 0]'),[0]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 1]'),[pi/2]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[-1 0]'),[pi]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 -1]'),[-pi/2]');
    %sph2cartm, type=0
    EXPECT_DOUBLE_EQ( [1 0]',sph2cartm([0 1]'));
    EXPECT_DOUBLE_EQ( [0 1]',sph2cartm([pi/2 1]'));
    EXPECT_DOUBLE_EQ( [-1 0]',sph2cartm([pi 1]'));
    EXPECT_DOUBLE_EQ( [0 -1]',sph2cartm([-pi/2 1]'));
    %sph2cartm2, type=0
    EXPECT_DOUBLE_EQ( [1 0]',sph2cartm2([0]'));
    EXPECT_DOUBLE_EQ( [0 1]',sph2cartm2([pi/2]'));
    EXPECT_DOUBLE_EQ( [-1 0]',sph2cartm2([pi]'));
    EXPECT_DOUBLE_EQ( [0 -1]',sph2cartm2([-pi/2]'));
    %sph2cartm/sph2cartm2, type=0
    val = randn(2,10); 
    EXPECT_NEAR( sph2cartm(cart2sphm(val)),val,1e-10 );
    val = randn(2,300); 
    val=normalizematrix(val,'colnorm',2);
    EXPECT_NEAR( sph2cartm2(cart2sphm2(val)), val, 1e-10 );

    %S^2
    %cart2sphm
    EXPECT_DOUBLE_EQ( cart2sphm([1 0 0]'),[0 0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 1 0]'),[pi/2 0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 0 1]'),[pi/2 pi/2 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([-1 0 0]'),[pi 0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 -1 0]'),[pi/2 pi 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 0 -1]'),[pi/2 -pi/2 1]');
    %type=2/3
    EXPECT_DOUBLE_EQ( cart2sphm([1 0 0]',2),[0 0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 1 0]',2),[pi/2 0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 0 1]',2),[0 pi/2 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 0 1]',3),[0 0 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([0 1 0]',3),[pi/2 pi/2 1]');
    EXPECT_DOUBLE_EQ( cart2sphm([1 0 0]',3),[0 pi/2 1]');
    %cart2sphm2, type=0
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[1 0 0]'),[0 0]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 1 0]'),[pi/2 0]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 0 1]'),[pi/2 pi/2]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[-1 0 0]'),[pi 0]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 -1 0]'),[pi/2 pi]');
    EXPECT_DOUBLE_EQ( cart2sphm2(2*[0 0 -1]'),[pi/2 -pi/2]');
    %sph2cartm, type=0
    EXPECT_DOUBLE_EQ( [1 0 0]',sph2cartm([0 0 1]'));
    EXPECT_DOUBLE_EQ( [0 1 0]',sph2cartm([pi/2 0 1]'));
    EXPECT_DOUBLE_EQ( [0 0 1]',sph2cartm([pi/2 pi/2 1]'));
    EXPECT_DOUBLE_EQ( [-1 0 0]',sph2cartm([pi 0 1]'));
    EXPECT_DOUBLE_EQ( [0 -1 0]',sph2cartm([pi/2 pi 1]'));
    EXPECT_DOUBLE_EQ( [0 0 -1]',sph2cartm([pi/2 -pi/2 1]'));
    %type=2/3
    EXPECT_DOUBLE_EQ( [1 0 0]',sph2cartm([0 0 1]',2) );
    EXPECT_DOUBLE_EQ( [0 1 0]',sph2cartm([pi/2 0 1]',2) );
    EXPECT_DOUBLE_EQ( [0 0 1]',sph2cartm([0 pi/2 1]',2) );
    EXPECT_DOUBLE_EQ( [0 0 1]',sph2cartm([0 0 1]',3) );
    EXPECT_DOUBLE_EQ( [0 0 1]',sph2cartm([pi/2 0 1]',3) );
    EXPECT_DOUBLE_EQ( [1 0 0]',sph2cartm([0 pi/2 1]',3) );
    %sph2cartm2, type=0
    EXPECT_DOUBLE_EQ( [1 0 0]',sph2cartm2([0 0]'));
    EXPECT_DOUBLE_EQ( [0 1 0]',sph2cartm2([pi/2 0]'));
    EXPECT_DOUBLE_EQ( [0 0 1]',sph2cartm2([pi/2 pi/2]'));
    EXPECT_DOUBLE_EQ( [-1 0 0]',sph2cartm2([pi 0]'));
    EXPECT_DOUBLE_EQ( [0 -1 0]',sph2cartm2([pi/2 pi]'));
    EXPECT_DOUBLE_EQ( [0 0 -1]',sph2cartm2([pi/2 -pi/2]'));
    %sph2cartm/sph2cartm2, type=0
    val = randn( 3, 300 );
    EXPECT_NEAR( sph2cartm(cart2sphm(val)), val, 1e-12);
    val = normalizematrix( val, 'colnorm',2 );
    EXPECT_NEAR( sph2cartm2(cart2sphm2(val)), val, 1e-12);
    %cart2sphm, ...type=1, test using Matlabs cart2sph()
    val = randn(3,300); 
    [val1,val2,val3] = cart2sph(val(1,:),val(2,:),val(3,:));
    EXPECT_DOUBLE_EQ( cart2sphm(val,2),[val1;val2;val3]);
    
    [val1,val2,val3]=sph2cart(val(1,:),val(2,:),val(3,:));
    EXPECT_DOUBLE_EQ( sph2cartm(val,2),[val1;val2;val3]);
    
    [val1,val2,~]=cart2sph(val(1,:),val(2,:),val(3,:));
    EXPECT_DOUBLE_EQ( cart2sphm2(val,2),[val1;val2]);
    
    [val1,val2,val3]=sph2cart(val(1,:),val(2,:),ones(1,size(val,2)));
    EXPECT_DOUBLE_EQ( sph2cartm2(val(1:2,:),2),[val1;val2;val3]);

    %S^99
    %sph2cartm/sph2cartm2, type=0
    val = randn(100,300); 
    EXPECT_NEAR( sph2cartm(cart2sphm(val)),val,1e-11);
    val = randn(100,300); 
    val=normalizematrix(val,'colnorm',2);
    EXPECT_NEAR( sph2cartm2(cart2sphm2(val)),val,1e-11);

    %sph2sphm
    val = sph2sphm( randn(100,300),0,0);   
    EXPECT_NEAR( val,sph2sphm(val,0,0),1e-9);
    val = sph2sphm2(randn(100,300),0,0);  
    EXPECT_NEAR( val,sph2sphm2(val,0,0),1e-10);
    %S^2, type 1
    val = sph2sphm( randn(3,300),2,2);   
    EXPECT_DOUBLE_EQ( val,sph2sphm(val,2,2) );
    val = sph2sphm2( randn(2,300),2,2);  
    EXPECT_DOUBLE_EQ( val,sph2sphm2(val,2,2) );
    val = sph2sphm( randn(3,300),0,0); 
    val1 = sph2sphm( sph2sphm(val,0,2),2,0);
    EXPECT_NEAR( val, val1, 1e-12 );
    val = sph2sphm( randn(3,300),2,2); 
    val1 = sph2sphm( sph2sphm(val,2,0),0,2);
    EXPECT_NEAR( val, val1, 1e-12 );
    %type 2
    val = sph2sphm( randn(3,300), 3, 3 );   
    EXPECT_NEAR( val,sph2sphm(val,3,3), 1e-10 );
    val = sph2sphm2(randn(2,300), 3, 3 );  
    EXPECT_NEAR( val,sph2sphm2(val,3,3), 1e-10 );
    val = sph2sphm( randn(3,300), 0, 0 ); 
    val1=sph2sphm( sph2sphm(val,0,2),2,0);
    EXPECT_NEAR( val, val1, 1e-10 );
    val = sph2sphm( randn(3,300), 3, 3 ); 
    val1=sph2sphm( sph2sphm(val,3,0), 0, 3 );
    EXPECT_NEAR( val, val1,1e-10 );
    %type 2 <-> type 3
    EXPECT_DOUBLE_EQ( sph2sphm([0 pi/2 1]',2,3), [0 0 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([2*pi pi/2 1]',2,3), [0 0 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([0 0 1]',2,3), [0 pi/2 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([0 0 1]',3,2), [0 pi/2 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([0 pi/2 1]',3,2), [0 0 1]' );
    %type 4 <-> type 2
    EXPECT_DOUBLE_EQ( sph2sphm([pi 0 1]',4,2), [-pi/2 0 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([pi/2 0 1]',4,2), [0 0 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([3*pi/2 0 1]',4,2), [-pi 0 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([0 pi/2 1]',4,2), [pi/2 -pi/2 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([pi/4 -pi/4 1]',4,2), [pi/4 pi/4 1]' );
    EXPECT_DOUBLE_EQ( sph2sphm([pi/4 pi/4 1]',4,2), [pi/4 -pi/4 1]' );

    %type 100 tests
    EXPECT_DOUBLE_EQ( cart2sphm([1;0;1],100), [0;1;1] );
    EXPECT_DOUBLE_EQ( [1;0;1], sph2cartm([0;1;1],100) );
    EXPECT_DOUBLE_EQ( cart2sphm2([2;0;2],100), [0;2] );
    EXPECT_DOUBLE_EQ( [1;0;1], sph2cartm2([0;1],100) );
    EXPECT_DOUBLE_EQ( cart2sphm([0;1;1],-100), [90;1;1] );
    EXPECT_DOUBLE_EQ( [0;1;1], sph2cartm([90;1;1],-100) );
    EXPECT_DOUBLE_EQ( cart2sphm2([0;2;2],-100), [90;2] );
    EXPECT_DOUBLE_EQ( [0;1;1], sph2cartm2([90;1],-100) );

    %dim 2 tests of type 2/3/4
    val = randn(2,100);
    EXPECT_NEAR( cart2sphm(val,0), cart2sphm(val,2), 1e-10 );
    EXPECT_NEAR( cart2sphm(val,0), cart2sphm(val,3), 1e-10 );
    EXPECT_DOUBLE_EQ( sph2cartm(val,0), sph2cartm(val,2), 10 );
    EXPECT_DOUBLE_EQ( sph2cartm(val,0), sph2cartm(val,3), 10 );
    EXPECT_DOUBLE_EQ( sph2cartm([0;1],4), [0;1], 10 );
    EXPECT_DOUBLE_EQ( cart2sphm([1;0],-2), [0;1], 10 );
    EXPECT_DOUBLE_EQ( cart2sphm([0;1],-2), [90;1], 10 );
    EXPECT_NEAR( cart2sphm(val,-1), cart2sphm(val,-2), 1e-6 );
    EXPECT_NEAR( cart2sphm(val,-1), cart2sphm(val,-3), 1e-6 );
    EXPECT_DOUBLE_EQ( sph2cartm(val,-1), sph2cartm(val,-2), 10  );
    EXPECT_DOUBLE_EQ( sph2cartm(val,-1), sph2cartm(val,-3), 10  );
    EXPECT_DOUBLE_EQ( sph2cartm([0;1],-4), [0;1], 10  );

    %tests of type 104
    EXPECT_DOUBLE_EQ( sph2cartm([pi/2;-1;20],104), [20;0;-1], 10 );
    EXPECT_DOUBLE_EQ( sph2cartm([90;-1;20],-104), [20;0;-1], 10 );

    %sph2sphm2 tests with TTEST_ALLTESTFLAG set missing


end

%% test convm
seq1 = convm( [1 2] );
if( INIT('all') )
    EXPECT_DOUBLE_EQ( seq1,[1 2]);
    seq1 = convm( [1 1],[1 2]);
    EXPECT_DOUBLE_EQ( seq1,[1 3 2]);
    seq1 = convm( [1 1]',[1 2]');
    EXPECT_DOUBLE_EQ( seq1,[1 3 2]');
    seq1 = convm( [1 2 1],[1; 3; 1],'mult',[2 1],'normalize');
    EXPECT_DOUBLE_EQ( seq1,1/80*[1 4 6 4 1;3 12 18 12 3;1 4 6 4 1]);
    seq1 = convm( [1 2 1],[1; 3; 1],'outer','normalize');
    EXPECT_DOUBLE_EQ( seq1,1/20*[1 3 1; 2 6 2; 1 3 1]);
    convm(randn(2,2,2),randn(2,2,2),randn(2,2,2),randn(2,2,2) );
    EXPECT_EQ( convm([],[1 2]), [] );
    EXPECT_EQ( convm([1 2],[]), [] );
    EXPECT_EQ( convm([1 2],sym([2 1])), sym([2 5 2]) );

    if(tmiscinstalled); 
        seq1 = convm(sym([1 1]),[1 2]); 
        EXPECT_DOUBLE_EQ( seq1,[1 3 2]);
        seq1 = convm(sym([1 1]),sym([1; 2])); 
        EXPECT_DOUBLE_EQ( seq1,[1 1;2 2]);
    end        
end

%% test dec2basem
num = dec2basem([10; 1],[1 2; -1 2]);
EXPECT_EQ( num, [1 1 2;1 1 1] );
num=dec2basem(-1,2,5);
EXPECT_EQ( num(1,1:5), [1 1 1 1 1] );
num=dec2basem(3,2,5);
EXPECT_EQ( num, [0 0 0 1 1] );
dec2basem( [10;-10],[2 1;0 2],5);
if( INIT('all') )
    EXPECT_ERROR( 'dec2basem(10,1/2,5);' );
end

%% test factorialm
EXPECT_EQ( factorialm([3 4]), 144 );

%% test gcdm
EXPECT_EQ( gcdm(10,15,20), 5 );
EXPECT_EQ( gcdm(1), 1);
EXPECT_EQ( gcdm(), 0);
EXPECT_EQ( gcdm([]), 0);
EXPECT_EQ( gcdm([10,15,20]), 5 );

%% test ind2subm');
ind2subm([2 1 3],[1 6]);
if( INIT('all') )
    EXPECT_ERROR( 'tryind2subm([2 1 3],[1 7]);' );
end



%% test isvectorm
isvectorval=[isvectorm(ones(1,1,3)) isvectorm(ones(0)) isvector(ones(1,1,3)) isvector(ones(0))];
EXPECT_EQ( isvectorval, [true false false false] );

%% test kronm
val=kronm([1 2 3],[2; 3]);
EXPECT_EQ( val, [2 4 6; 3 6 9] );
val=kronm([],[],[]);
EXPECT_EQ( val, [] );

%% test lcmm
EXPECT_EQ( lcmm(10,15,20), 60 );
EXPECT_EQ( lcmm([10,15,20]), 60 );
EXPECT_EQ( lcmm(1), 1 );
EXPECT_EQ( lcmm(), 1 );
EXPECT_EQ( lcmm(-4), 4 );
EXPECT_EQ( lcmm(-4,0), 0 );
EXPECT_EQ( lcmm(0), 0 );
EXPECT_EQ( lcmm([]), 1 );


%% test maxm
EXPECT_EQ( maxm(), -Inf );
EXPECT_EQ( maxm([1 2 2; 2 3 4; 2 1 3;4 5 NaN]), 5);
EXPECT_EQ( maxm([1 2 2; 2 3 4; 2 1 3],1), [2 3 4] );
EXPECT_EQ( maxm([1 0; 1 -1],[]), [1 0;1 1] );

%% test minm
EXPECT_EQ( minm(), Inf );
EXPECT_EQ( minm([1 2 2; 2 3 4; 2 1 3;4 5 NaN]), 1 );
EXPECT_EQ( minm([1 2 2; 2 3 4],2), [1;2] );
EXPECT_EQ( minm([1 0; 1 -1],[]), [1 0;1 1] );

%% test nchoosekm
EXPECT_EQ( nchoosekm([10 5],[5 2]),2520 );

%% test ndimsm
EXPECT_EQ( ndimsm([1; 2; 3]),1);

%% test onesm
EXPECT_EQ( onesm,[]);
EXPECT_EQ( onesm([]),[]);
EXPECT_EQ( onesm(1,2),[1 1]);
EXPECT_EQ( onesm(3),[1;1;1]);

%% test padarraym
EXPECT_EQ( padarraym([1],1),[0;1;0]);
if( INIT('all') )
    EXPECT_EQ( padarraym([],2), [] );
    EXPECT_EQ( padarraym([1;2],1), [0;1;2;0] );
    EXPECT_EQ( padarraym({1;2},1), {[];1;2;[]} );
    EXPECT_EQ( padarraym([1 2],1,'pre'), [0 0 0;0 1 2] );
    EXPECT_EQ( padarraym({1 2},1,'pre'), {[] [] [];[] 1 2} );
    EXPECT_EQ( padarraym([1 2],1,'post'), [1 2 0;0 0 0] );
    EXPECT_EQ( padarraym({1 2},1,'post'), {1 2 [];[] [] []} );
    EXPECT_EQ( padarraym([1 2],1), [0 0 0 0;0 1 2 0;0 0 0 0] );
    EXPECT_EQ( padarraym({1 2},1), {[] [] [] [];[] 1 2 [];[] [] [] []} );
    EXPECT_EQ( padarraym([1 2;3 4],0), [1 2;3 4] );
    EXPECT_EQ( padarraym({1 2;3 4},0), {1 2;3 4} );
    EXPECT_EQ( padarraym({1 2;@sum 4},1), {[] [] [] [];[] 1 2 [];[] @sum 4 [];[] [] [] []});
    EXPECT_ERROR( 'padarraym([1;2],-1);' );
end

%% test parsem
[v,arg]=parsem('selftest',{'asd','selftest',2},@setupm);
EXPECT_EQ( v, 2 );
EXPECT_EQ(arg, {'asd'});
if( INIT('all') )
    [v,arg]=parsem('selftest','selftest');
    EXPECT_EQ( v,1);

    [v,arg]=parsem('tommsch',{'asd','tommsch',2},'is the best','set');
    EXPECT_FALSE( ~isequal( v,'is the best') || ~isequal(arg, {'asd', 'tommsch', 'is the best'}));

    [v,arg]=parsem('tommsch',{'asd','tommsch',0,'tommsch',1},'is the second best','condset');
    EXPECT_FALSE( ~isequal( v,1) || ~isequal(arg, {'asd', 'tommsch', 1}));

    [v,arg]=parsem('tommsch',{'asd'});
    EXPECT_FALSE( ~isequal( v,0) || ~isequal(arg, {'asd'})); 

    [v,arg]=parsem({'tommsch','tomm'},{'tomm',10,'asd','tommsch',2},'is the best','set');
    EXPECT_FALSE( ~isequal( v,'is the best') || ~isequal(arg, {'asd', 'tommsch', 'is the best'}));

    [v,arg]=parsem({'tomm','tommsch'},{'asd','tommsch',0,'tommsch',1},'is the second best','condset');
    EXPECT_FALSE( ~isequal( v,1) || ~isequal(arg, {'asd', 'tomm', 1}));

    [v,arg]=parsem({'t','tommsch','t'},{'asd'});
    EXPECT_FALSE( ~isequal( v,0) || ~isequal(arg, {'asd'}));

    [v,arg]=parsem({},{'asd'});
    EXPECT_FALSE( ~isequal( v,0) || ~isequal(arg, {'asd'}));

    [v,arg]=parsem({},{});
    EXPECT_FALSE( ~isequal( v,0) || ~isequal(arg, {}));

    EXPECT_WARNING( 'parsem(''asd'',{''asd''},1);', 'parsem:missing' );
    [v,arg]=parsem({},{'parse_one'});
    EXPECT_EQ( v,1);
    [v,arg]=parsem({},{'parse_random'});
    EXPECT_TRUE( isequal(v,0) || isequal(v,1) );
    parsem('...',{'help'}); 
    fprintf( '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b' ); %to remove line breaks and other white space printed by parsem('help')
    EXPECT_WARNING( 'parsem({''2 is the oddest prime.''},''test'');', 'parsem:unkown' ); %make warnings to error

    EXPECT_WARNING( 'parsem({''name''},''name'',''val'');', 'parsem:missing' ); %make warnings to error        

    [v,arg,ret] = parsem('st',{'selftest','st'});
    EXPECT_EQ( ret,'st');
    [v,arg,ret] = parsem({'st','selftest',},{'selftest'});
    EXPECT_EQ( ret,'selftest');
    [v,arg,ret] = parsem({'st','selftest',},{'t'});
    EXPECT_EQ( ret,'st');

    warning( 'error', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',11}, ''expect'', {''clop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',10}, ''expect'', {''clop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',-1}, ''expect'', {''clop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',[]}, ''expect'', {''clop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',11}, ''expecte'',{''clop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',10}, ''expecte'',{''clop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',-1}, ''expecte'',{''clop'',[0 10]} );', 'parsem:expect' );
    parsem( 't', {'t',[]}, 'expecte',{'clop',[0 10]} );

    EXPECT_WARNING( 'parsem( ''t'', {''t'',11}, ''expect'',{''clcl'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',-1}, ''expect'',{''clcl'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',[]}, ''expect'',{''clcl'',[0 10]} );', 'parsem:expect' );
    parsem( 't', {'t',10}, 'expect',{'clcl',[0 10]} );        

    EXPECT_WARNING( 'parsem( ''t'', {''t'',11}, ''expect'',{''opcl'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',0},  ''expect'',{''opcl'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',-1}, ''expect'',{''opcl'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',[]}, ''expect'',{''opcl'',[0 10]} );', 'parsem:expect' );
    parsem( 't', {'t',9}, 'expect',{'opcl',[0 10]} );     

    EXPECT_WARNING( 'parsem( ''t'', {''t'',11}, ''expect'',{''opop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',10}, ''expect'',{''opop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',0},  ''expect'',{''opop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',-1}, ''expect'',{''opop'',[0 10]} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',[]}, ''expect'',{''opop'',[0 10]} );', 'parsem:expect' );
    parsem( 't', {'t',9}, 'expect',{'opcl',[0 10]} );          

    EXPECT_WARNING( 'parsem( ''t'', {''t'',11}, ''expect'',{0 1 10} );', 'parsem:expect' );
    EXPECT_WARNING( 'parsem( ''t'', {''t'',[]}, ''expect'',{0 1 10} );', 'parsem:expect' );
    parsem( 't', {'t',[]}, 'expecte',{0 1 10} );
    parsem( 't', {'t',10}, 'expect',{0 1 10} );

    EXPECT_WARNING( 'parsem( ''t'', {''t'',@sum}, ''expect'',@isnumeric );', 'parsem:expect' );
    parsem( 't', {'t',10}, 'expect',@isnumeric );

    EXPECT_WARNING( 'parsem( {''s'',''t'',''u''}, {''v''}, ''expect'',''name'' );', 'parsem:expect' );
    parsem( {'s','t','u'}, {'t'}, 'expect','name' );    

end

%% test plotm
try; %#ok<TRYNC>
    evalc('plot([2]);'); %trigger possible warning: Warning: MATLAB has disabled some advanced graphics rendering features by switching to software OpenGL....
    close all; end;

EXPECT_NO_THROW( 'plotm(randn(1,40));' );
if( INIT('all')  );
    %1d-stuff


    EXPECT_WARNING( 'plotm(randn(1,40),''wrongoption'');', 'plotm:option' );
    EXPECT_NO_THROW( 'plotm(randn(2,40),''.-'',''box'',2); title(''plotm: Some lines and something like a square'');' );
    EXPECT_WARNING( 'plotm(num2cell(reshape(cumsum(rand(20,1)),2,[])'',2),''wrongoption'',1);', 'plotm:option' );
    EXPECT_NO_THROW( 'plotm(num2cell(reshape((rand(4,1)),2,[])'',2),''verbose'',1);' );
    EXPECT_NO_THROW( 'plotm(num2cell(reshape(cumsum(rand(20,1)),2,[])'',2),''verbose'',1);' );
    if(sequenceinstalled);
        EXPECT_NO_THROW( 'plotm(sequence(randn(20,1),[0])); ' );
    end
    %2d-stuff      

    EXPECT_NO_THROW( 'plotm(randn(1,10)+1i*randn(1,10),''.'');');

    val = randn(2,20); 
    EXPECT_NO_THROW( 'plotm(val,''.''); ');
    EXPECT_NO_THROW( 'plotm(val,''boundary'',.5,''Color'',''red'');');
    EXPECT_NO_THROW( 'plotm(val,''boundary'',-.5,''Color'',''blue''); ');
    EXPECT_NO_THROW( 'plotm(val,''boundary'',-1.5,''Color'',''green'');');
    EXPECT_NO_THROW( 'plotm(val,''hull'',''Color'',''black'');');

    if(sequenceinstalled);
        EXPECT_NO_THROW( 'plotm(sequence([1 1; 1 0;3 1],[0 0]));');
    end

    %3d-stuff
    EXPECT_NO_THROW( 'plotm( randn(3,40), ''resolution'',0, ''MarkerSize'',100 ); ');
    EXPECT_NO_THROW( 'title(''plotm: A pointcloud and a box.'');');
    EXPECT_NO_THROW( 'plotm(''box'',3);');

    EXPECT_NO_THROW( 'plotm(randn(3,40),''resolution'',100,''surface'',''contour''); view(0,50);');
    EXPECT_NO_THROW( 'plotm([0 0 0;1 0 0;1 1 0;0 1 0; 0 0 0;0 0 1;1 0 1;1 1 1;0 1 1; 0 0 1;0 0 0; 1 0 0;1 0 1;0 0 1; 0 0 0;0 1 0;1 1 0;1 1 1;0 1 1; 0 1 0;0 0 0;0 1 0;0 1 1;0 0 1;0  0 0;1 0 0;1 1 0;1 1 1;1 0 1;1  0 0;0 0 0]'',''.-'');');
    if(sequenceinstalled && tmiscinstalled)
        EXPECT_NO_THROW( 'plotm(sequence(randn(20),[0;0]));');
    end

    val = randn(3,20);
    EXPECT_NO_THROW( 'plotm(val,''boundary'',.5);');
    EXPECT_NO_THROW( 'plotm(val,''hull'');');

    %4d-stuff
    val = randn(4,20);
    EXPECT_NO_THROW( 'plotm(val);');

    %400d-stuff
    val = randn(400,2);
    EXPECT_WARNING( 'plotm(val);', 'plotm:dim' );
end;

%% test repcellm
EXPECT_EQ( repcellm(10 ,2),{10;10});
EXPECT_EQ( repcellm(10 ,[2 3]),{10 10 10;10 10 10});

%% test repmatm
EXPECT_EQ( repmatm(10 ,2), [10;10] );
EXPECT_EQ( repmatm(10 ,[2 3]), [10 10 10;10 10 10] );
EXPECT_EQ( repmatm([1 2] ,2, 3), [1 2 1 2 1 2;1 2 1 2 1 2] );

%% test sizem
[a,b,c]=sizem(zeros(2,0,1));
EXPECT_EQ( a,2);
EXPECT_EQ(b,0);
EXPECT_EQ(c,1);
if( INIT('all') )
    a=sizem(1); EXPECT_EQ( a,1);
    a=sizem([1;2]); EXPECT_EQ( a,2);
    a=sizem([1 2]); EXPECT_EQ( a,[1 2]);
    a=sizem([]); EXPECT_EQ( a,[]);
    EXPECT_EQ( 3,sizem([1 2;3 4;5 6],1));
end

%% test sph2cartm % tested above
%% test sph2cart2m % tested above
%% test sph2sph % tested above
%% test sph2sphm % tested above

%% test squeezem
EXPECT_EQ( squeezem([2 1 3]),[2;1;3]);

%% test summ
EXPECT_EQ( summ([2 3; -4 -5],[],'abs'),14);
EXPECT_EQ( summ([2 3; -4 -5]),-4);
EXPECT_EQ( summ([2 3; -4 -5],1),[-2 -2]);

%% test upsamplem
try
    [d,dmin]=upsamplem([1 2 ; 3 4], [0;0], [2 1; 0 -2],Inf);
    EXPECT_EQ( d,[inf inf 1;2 inf inf;inf inf 3;4 inf inf], dmin,[0;-2] );

    [d,dmin]=upsamplem([1 2 ; 3 4], [], [2 1; 0 -2],Inf);
    EXPECT_EQ( d,[inf inf 1;2 inf inf;inf inf 3;4 inf inf], dmin,[0;-2] );
    if(INIT('all') )
        [d,dmin]=upsamplem([],0,2);
        EXPECT_EQ( d,[], dmin,[0]);
        [d,dmin]=upsamplem([1 2 3]',0,-1);
        EXPECT_EQ( d,[3 2 1]', dmin,-2);
        [d,dmin]=upsamplem([1 2 3]',0,1);
        EXPECT_EQ( d,[1 2 3]', dmin,0);
        [d,dmin]=upsamplem([1 2 ; 3 4], [0;0], [2 0;0 2]);
        EXPECT_EQ( d,[1 0 2;0 0 0;3 0 4], dmin,[0;0]);
        
        EXPECT_ERROR( 'upsamplem([1 2 3],[0;0]);');
        EXPECT_ERROR( 'upsamplem([1 2 3],[0;0],[2]);');
        EXPECT_ERROR( 'upsamplem([1 2 3],[2],[1 2]);');
        EXPECT_ERROR( 'upsamplem([1 2 3],[2 2],[2 0;0 2]);');
        EXPECT_ERROR( 'upsamplem([1 2 3],[2;2],[2 0;0 2],[1 0]);');
    end
catch
    fprintf( 'ERROR: ''upsamplem'' broken. Probably tmisc package not installed.\n');
end

%% test zerosm
EXPECT_EQ( zerosm(),[]);
EXPECT_EQ( zerosm([]),[]);
EXPECT_EQ( zerosm([1,2]),[0 0]);
EXPECT_EQ( zerosm(3),[0;0;0]);
EXPECT_EQ( zerosm(1,1,2),permute([0; 0],[3 2 1]));
