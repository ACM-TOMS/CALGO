% test TTEST

INIT('clear');
INIT( 'exp',@()(0) );
INIT( 'ass',@()(0) );
INIT( 'v',0 );

% preconditions   

%% test FAIL

assert( ~EXPECT_FAIL );

%% test SUCCEED
assert( EXPECT_SUCCEED );
assert( ASSERT_SUCCEED );

%% test TRUE
assert( EXPECT_TRUE(1) );
assert( EXPECT_TRUE(2) );
assert( ~EXPECT_TRUE(0) );
assert( ASSERT_TRUE(1) );
assert( ASSERT_TRUE(2) );
assert( ~ASSERT_TRUE(0) );


%% test FALSE
assert( EXPECT_FALSE(0) );
assert( ~EXPECT_FALSE(1) );
assert( ASSERT_FALSE(0) );
assert( ~ASSERT_FALSE(1) );

%% test EQ
assert( EXPECT_EQ([],[]) );
assert( EXPECT_EQ(2,2) );
assert( EXPECT_EQ(@sum,@sum) );
assert( ~EXPECT_EQ(@sum,@prod) );
assert( ~EXPECT_EQ(2,3) );

%% test NE
assert( ~EXPECT_NE(2,2) );
assert( ~EXPECT_NE(@sum,@sum) );
assert( EXPECT_NE(@sum,@prod) );
assert( EXPECT_NE(2,3) );

%% test LE
assert( EXPECT_LE(2,2) );
assert( EXPECT_LE(2,3) );
assert( ~EXPECT_LE(3,2) );

%% test LT
assert( ~EXPECT_LT(2,2) );
assert( EXPECT_LT(2,3) );
assert( ~EXPECT_LT(3,2) );

%% test GE
assert( EXPECT_GE(2,2) );
assert( ~EXPECT_GE(2,3) );
assert( EXPECT_GE(3,2) );

%% test GT
assert( ~EXPECT_GT(2,2) );
assert( ~EXPECT_GT(2,3) );
assert( EXPECT_GT(3,2) );

%% test STREQ
assert( EXPECT_STREQ('23','23') );
assert( EXPECT_STREQ('23','23') );
assert( EXPECT_STREQ('ab','ab') );
assert( ~EXPECT_STREQ('ab','aB') );
assert( ~EXPECT_STREQ('ab','ac') );
assert( ~EXPECT_STREQ('23','') );
assert( EXPECT_STREQ('','') );

%% test STRNEQ
assert( ~EXPECT_STRNE('23','23') );
assert( ~EXPECT_STRNE('23','23') );
assert( ~EXPECT_STRNE('ab','ab') );
assert( EXPECT_STRNE('ab','aB') );
assert( EXPECT_STRNE('ab','ac') );
assert( EXPECT_STRNE('23','') );
assert( ~EXPECT_STRNE('','') );

%% test STRCASEEQ
assert( EXPECT_STRCASEEQ('23','23') );
assert( EXPECT_STRCASEEQ('23','23') );
assert( EXPECT_STRCASEEQ('ab','ab') );
assert( EXPECT_STRCASEEQ('ab','aB') );
assert( ~EXPECT_STRCASEEQ('ab','ac') );
assert( ~EXPECT_STRCASEEQ('23','') );
assert( EXPECT_STRCASEEQ('','') );

%% test STRCASENEQ
assert( ~EXPECT_STRCASENE('23','23') );
assert( ~EXPECT_STRCASENE('23','23') );
assert( ~EXPECT_STRCASENE('ab','ab') );
assert( ~EXPECT_STRCASENE('ab','aB') );
assert( EXPECT_STRCASENE('ab','ac') );
assert( EXPECT_STRCASENE('23','') );
assert( ~EXPECT_STRCASENE('','') );

%% test PRED
assert( EXPECT_PRED( @isempty, [] ) );
assert( EXPECT_PRED( @isempty, {} ) );
assert( ~EXPECT_PRED( @isempty, 1 ) );
assert( EXPECT_PRED( @(a,b) a<b, 1,2 ) );
assert( EXPECT_PRED( @(a,b,c) a+b+c==0, 1,2,-3 ) );

%% test NPRED
assert( ~EXPECT_NPRED( @isempty, [] ) );
assert( ~EXPECT_NPRED( @isempty, {} ) );
assert( EXPECT_NPRED( @isempty, 1 ) );
assert( ~EXPECT_NPRED( @(a,b) a<b, 1,2 ) );
assert( ~EXPECT_NPRED( @(a,b,c) a+b+c==0, 1,2,-3 ) );

%% test NEAR
assert( EXPECT_NEAR( 1, 2, 1.1 ) );
assert( ~EXPECT_NEAR( 1, 2, 0.9 ) );

%% test DOUBLE_EQ
assert( EXPECT_DOUBLE_EQ( 1, 1+eps ) );
assert( EXPECT_DOUBLE_EQ( 1, 1+4*eps, 4 ) );
assert( ~EXPECT_DOUBLE_EQ( 1, 1+5*eps, 4 ) );

%% test SINGLE_EQ
assert( EXPECT_SINGLE_EQ( single(1), 1+eps(single(1)) ) );
assert( ~EXPECT_SINGLE_EQ( single(1), 1+5*eps(single(1)), 4 ) );

%% test NO_THROW
assert( EXPECT_NO_THROW('2;') );
assert( ~EXPECT_NO_THROW('adsasfasf a asdkasd;') );
assert( EXPECT_NO_THROW('adsasfasf a asdkasd;','asd:asd', 'MATLAB:UndefinedFunction') );

%% test WARNING
assert( ~EXPECT_WARNING('0;') );
assert( EXPECT_WARNING('inv(0);') );
assert( ~EXPECT_WARNING('inv(0);', 'asd:asd') );
assert( EXPECT_WARNING('inv(0);', 'asd:asd', 'MATLAB:singularMatrix') );

%% test ERROR
assert( ~EXPECT_ERROR('0;') );
assert( EXPECT_ERROR('adsasfasdasda ajkdwaklsd;') );
assert( ~EXPECT_ERROR('adsasfasdasda ajkdwaklsd;', 'asd:asd' ) );
assert( EXPECT_ERROR('adsasfasdasda ajkdwaklsd;', 'asd:asd', 'MATLAB:UndefinedFunction' ) );

%% test R2017b
try;
    val = runtests('testTTESTR2017b');
    if( ~all([val.Passed]) );
        assert( false ); end;
catch
    assert( false );
end;
