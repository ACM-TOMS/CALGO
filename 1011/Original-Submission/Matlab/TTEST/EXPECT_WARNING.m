%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = EXPECT_WARNING( varargin )

global TTEST_ERROR_EXPECT;
global TTEST_ERRORFLAG;

ret = true;
try;
    cmd = varargin{1};
    s = warning; %save warning settings
    lastwarn( 'TTEST_NOWARN' );
    for i = 2:nargin; %make warnings to error
        warning( 'error', varargin{i} ); end; %#ok<CTPCT>

    if( nargin>1 );
        id = varargin(2:end);
    else;
        id = {}; end;
    %check input
    if( ~isstring(cmd) && ~ischar(cmd) );
        error( 'TTEST_EVAL:format', '\ncmd must be a string passable to evalin()' ); end;

    try;
        evalc( 'evalin( ''caller'', cmd );' ); %execute command
    if( nargin>1 );
        TTEST_FPRINTF( '\nExpected warning not thrown:\nExpected Warning: ' );
        TTEST_DISP( varargin(2:end) );
        if( ~isequal(lastwarn,'TTEST_NOWARN') );
            [a,b] = lastwarn;
            TTEST_FPRINTF( '\nThrown warning: %s:%s\n', b,a); end;
        ret = false; end;
    catch me;
        lastwarn( 'TTEST_ERROR' );
        if( ~isempty(id) );
            switch me.identifier
                case id; %check if at least one correct warning is thrown
                    %do nothing
                otherwise;
                    ret = false;
                    TTEST_FPRINTF( '\nFalse warning thrown:\n   Expected Warning: ' );
                    TTEST_DISP( varargin(2:end) );
                    TTEST_FPRINTF( '   Thrown warning: ');
                    TTEST_DISP( me.identifier ); end;
        else;
            TTEST_FPRINTF( '\nFalse warning/error thrown: %s:%s\n', me.identifier, me.message );
            ret = false; end; end;

    if( isequal(lastwarn,'TTEST_NOWARN') );
        TTEST_FPRINTF( '\n Expected warning, but no warning was thrown.\n ' );
        ret = false; end;

    warning( s ); % redo warning settings
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

