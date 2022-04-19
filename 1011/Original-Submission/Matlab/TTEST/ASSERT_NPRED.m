%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = ASSERT_NPRED( varargin )

global TTEST_ERROR_ASSERT;
global TTEST_ERRORFLAG;

ret = true;
try;
    pred = varargin{1};
    val = varargin(2:end);
    if( pred(val{:}) );
        ret = false;
        TTEST_FPRINTF( 'Value(s) do(es) fulfil predicate (but it/they should), where' );
        TTEST_FPRINTF( '\nPredicate = ' );
        TTEST_DISP( pred );
        for ii = 1:numel( val );
            TTEST_FPRINTF( ['\nValue ' num2str(ii) ' = '] );
            TTEST_DISP( val{ii} ); end; end;
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

