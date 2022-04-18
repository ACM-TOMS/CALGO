% test double

%#ok<*CTPCT>
%#ok<*ASGLU>
%#ok<*NASGU>

% preconditions    

%% test isAlways
EXPECT_EQ( 1,isAlways(4));
EXPECT_EQ( 0,isAlways(0));
EXPECT_EQ( 1,isAlways(-1));
EXPECT_EQ( 0,isAlways(NaN,[],'false'));
EXPECT_EQ( 1,isAlways(NaN,[],'true'));
EXPECT_EQ( 1,isAlways(Inf)); 
if( INIT('all') )
    EXPECT_WARNING( 'isAlways(NaN)','double:isAlways' );
    EXPECT_EQ(isAlways(NaN,[],'false'),0);
end

%% test simplify
EXPECT_EQ( 4,simplify(4));
EXPECT_PRED( @isnan, simplify(NaN));
EXPECT_EQ( inf,simplify(inf));
EXPECT_EQ( 1/0.,simplify(1/0.));
