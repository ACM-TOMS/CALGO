%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = ASSERT_ERROR( varargin )

global TTEST_ERROR_ASSERT;
global TTEST_ERRORFLAG;

ret = true;
try;
    cmd = varargin{1};
    %check input
    if( ~isstring(cmd) && ~ischar(cmd) );
        error( 'TTEST_EVAL:format', '\ncmd must be a string passable to evalin()' ); end;
    if( nargin>1 );
        id = varargin(2:end);
    else;
        id = {}; end;
    try;
        evalin( 'caller', cmd );
        TTEST_FPRINTF( '\nNo Error was thrown, but error was expected.\n');
        ret = false;
    catch me;
        if( ~isempty(id) );
            id = varargin(2:end);
            switch me.identifier
                case id;
                    % do nothing
                otherwise;
                    ret = false;
                    TTEST_FPRINTF( '\nFalse error thrown.\n   Thrown Error: %s\n   Expected error: %s\n', me.identifier, cmd ); end; end; end;
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

