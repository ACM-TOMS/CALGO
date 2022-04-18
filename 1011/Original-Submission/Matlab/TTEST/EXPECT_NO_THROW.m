%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = EXPECT_NO_THROW( varargin )

global TTEST_ERROR_EXPECT;
global TTEST_ERRORFLAG;

ret = true;
try;
    cmd = varargin{1};
    if( nargin>1 );
        id = varargin(2:end);
    else;
        id = {}; end;
    lastwarn( 'TTEST_NOWARN' );
    if( ~isstring(cmd) && ~ischar(cmd) );
        error( 'TTEST_EVAL:format', '\ncmd must be a string passable to evalin()' ); end;

    try;
        evalin( 'caller', cmd );
    catch Me;
        %error was thrown
        switch Me.identifier
            case id;
                %do nothing
            otherwise;
                ret = false;
                TTEST_FPRINTF( '\nError was thrown, but no error was expected. \n  Error: %s\n', Me.identifier ); end; end;

    %check if warning was thrown
    [warnMsg,warnID] = lastwarn;
    switch warnID;
        case id;
            %do nothing
        otherwise;
            if( ~isequal(warnMsg,'TTEST_NOWARN') );
                ret = false;
                TTEST_FPRINTF( '\nWarning was thrown, but no warning was expected. \n  Warning: %s\n  Warning ID: %s\n', warnMsg, warnID ); end; end;
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

