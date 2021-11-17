% test TTEST R2017b

% Tests which need at least R2017b to work
% this function must be called from testTTEST
% since there is no INIT part

% preconditions   

%% test EQ
assert( EXPECT_EQ('a','a') );
assert( EXPECT_EQ("a","a") );
assert( ~EXPECT_EQ('a',"a") );
assert( ~EXPECT_EQ("a",'a') );

%% test NE
assert( ~EXPECT_NE('a','a') );
assert( ~EXPECT_NE("a","a") );
assert( EXPECT_NE('a',"a") );
assert( EXPECT_NE("a",'a') );

%% test STREQ
assert( EXPECT_STREQ('23',"23") );
assert( EXPECT_STREQ('ab',"ab") );
assert( ~EXPECT_STREQ('ab',"aB") );
assert( ~EXPECT_STREQ('ab',"ac") );
assert( EXPECT_STREQ("23","23") );
assert( EXPECT_STREQ("ab","ab") );
assert( ~EXPECT_STREQ("ab","aB") );
assert( ~EXPECT_STREQ("ab","ac") );


%% test STRNEQ
assert( ~EXPECT_STRNE('23',"23") );
assert( ~EXPECT_STRNE('ab',"ab") );
assert( EXPECT_STRNE('ab',"aB") );
assert( EXPECT_STRNE('ab',"ac") );
assert( ~EXPECT_STRNE("23","23") );
assert( ~EXPECT_STRNE("ab","ab") );
assert( EXPECT_STRNE("ab","aB") );
assert( EXPECT_STRNE("ab","ac") );

%% test STRCASEEQ
assert( EXPECT_STRCASEEQ('23',"23") );
assert( EXPECT_STRCASEEQ('ab',"ab") );
assert( EXPECT_STRCASEEQ('ab',"aB") );
assert( ~EXPECT_STRCASEEQ('ab',"ac") );
assert( EXPECT_STRCASEEQ("23","23") );
assert( EXPECT_STRCASEEQ("AB","ab") );
assert( EXPECT_STRCASEEQ("ab","aB") );
assert( ~EXPECT_STRCASEEQ("ab","ac") );

%% test STRCASENEQ
assert( ~EXPECT_STRCASENE('23',"23") );
assert( ~EXPECT_STRCASENE('ab',"ab") );
assert( ~EXPECT_STRCASENE('ab',"aB") );
assert( EXPECT_STRCASENE('ab',"ac") );
assert( ~EXPECT_STRCASENE("23","23") );
assert( ~EXPECT_STRCASENE("AB","ab") );
assert( ~EXPECT_STRCASENE("ab","aB") );
assert( EXPECT_STRCASENE("ab","ac") );

%% test NO_THROW
assert( EXPECT_NO_THROW("2;") );
assert( ~EXPECT_NO_THROW("adsasfasf a asdkasd;") );
assert( ~EXPECT_NO_THROW("adsasfasf a asdkasd;","asd:asd", "MATLAB:UndefinedFunctio") );

%% test WARNING
assert( ~EXPECT_WARNING("0;") );
assert( EXPECT_WARNING("inv(0);") );
assert( ~EXPECT_WARNING("inv(0);", "asd:asd") );
assert( EXPECT_WARNING("inv(0);", "asd:asd", "MATLAB:singularMatrix") );

%% test ERROR
assert( ~EXPECT_ERROR("0;") );
assert( EXPECT_ERROR("adsasfasdasda ajkdwaklsd;") );
assert( ~EXPECT_ERROR("adsasfasdasda ajkdwaklsd;", "asd:asd" ) );
assert( EXPECT_ERROR("adsasfasdasda ajkdwaklsd;", "asd:asd", "MATLAB:UndefinedFunction" ) );
