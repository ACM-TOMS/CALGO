% testfile for example_ttest()

testflag = INIT( 'testflag' ); %Read out variable testflag

% preconditions    


%% test two
EXPECT_EQ(2,3);
EXPECT_EQ(2,1+1);

%% test three
EXPECT_DOUBLE_EQ( [3 4;1 2], [3+eps(1000) 4; 1 2] );
if( testflag ); %if testflag is set, make more tests
    EXPECT_DOUBLE_EQ( 3, 3+eps(3) );
    EXPECT_GT( 3, 3 );
end;

ASSERT_EQ(3,2); %this test triggers a matlab-assert, thus
EXPECT_EQ(3,2); %this statement is not executed anymore
    
%% test zero
EXPECT_PRED( @isscalar, 0 );
EXPECT_WARNING( "inv(0);", "wrong:identifier" );
b = [0 0];
a = [0];
EXPECT_ERROR( "[a;b];" );

