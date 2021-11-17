%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = EXPECT_FALSE( varargin )

global TTEST_ERROR_EXPECT;
global TTEST_ERRORFLAG;

ret = true;
try;
   if( varargin{1} );
       TTEST_FPRINTF( '\nArgument is true (but should be false), where' );
       TTEST_FPRINTF( '\nargument = ' );
       TTEST_DISP( varargin{1} );
       ret = false; end;
catch;
    error( 'TTEST_EVAL:format', 'Wrong format of varargin' ); end;

if( isequal(ret,false) );
    TTEST_ERRORFLAG = true;

    try;
        if( isequal(TTEST_ERROR_EXPECT,'exp') );
            TTEST_ERROR_EXPECT_WORKER();
        elseif( isequal(TTEST_ERROR_EXPECT,'ass') );
            TTEST_ERROR_ASSERT_WORKER();
        else;
            feval( TTEST_ERROR_EXPECT ); end;
    catch me;
        if( isequal(me.identifier,'MATLAB:assertion:failed') );
            assert( false ); 
        else; 
            rethrow( me ); end; end; end;
end

