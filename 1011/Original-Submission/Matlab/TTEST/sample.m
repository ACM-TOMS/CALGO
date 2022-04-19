function sample( testflag )
% This is an example file showing how to use the TTESTs

clc

%parse input
if( nargin<=0 || isempty(testflag) );
    testflag = 0; end;

%pre-processing
INIT( 'clear' );               % Reset all global variables from TTEST SUITE
INIT( 'testflag',testflag );   % Set persistent variable testflag in INIT()
                               % This is a way to pass variables to a script test suite

%start test
runtests( 'testsample' );

%parse return values of test
if( INIT('errorflag') );                       % Check TTEST_ERRORFLAG, which is set if any of the tests fails.
    fprintf( 'Test failed\n' ); end;           %   EXPECT_XXX statements do not trigger a matlab-assert(), 
                                               %   thus, if such tests fail it is not recognized by Matlab.
                                               %   and so we test it by hand
                                               

