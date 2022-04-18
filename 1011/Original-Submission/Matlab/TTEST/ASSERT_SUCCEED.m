%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = ASSERT_SUCCEED( varargin )

global TTEST_ERROR_ASSERT;
global TTEST_ERRORFLAG;

ret = true;
try;
    ret = true;
catch;
    error( 'TTEST_EVAL:format', 'Wrong format of varargin' ); end;

if( isequal(ret,false) );
    TTEST_ERRORFLAG = true;

    try;
        if( isequal(TTEST_ERROR_ASSERT,'exp') );
            TTEST_ERROR_EXPECT_WORKER();
        elseif( isequal(TTEST_ERROR_ASSERT,'ass') );
            TTEST_ERROR_ASSERT_WORKER();
        else;
            feval( TTEST_ERROR_ASSERT ); end;
    catch me;
        if( isequal(me.identifier,'MATLAB:assertion:failed') );
            assert( false ); 
        else; 
            rethrow( me ); end; end; end;
end

