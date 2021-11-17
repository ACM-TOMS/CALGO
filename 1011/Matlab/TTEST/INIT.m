function ret = INIT( varargin )
% INIT
% ret = INIT( 'name',[val] )
% ret = INIT( [option] );
% The main purpose of this funciton is to set the environment for a new test session.
%
% INIT( 'name' )
% INIT( 'name',val )
% If 'name' is not any of the following: 'flag','clear','v','verbose','p','print' then
%   ret = INIT( 'name' )
%       Creates a new persistant variable with name 'name' and value [] OR
%       returns the value of the persistant variable with name 'name'
%
%   ret = INIT( 'name',val )
%       Creates a new persistant variable with name 'name' and value 'val' OR
%       returns the saved value of the persistant variable with name 'name' and sets the new value to 'val'
%
% ret = INIT( 'flag',[val] )
%   Sets the value of the global variable TTEST_ERRORFLAG to val (if val is given) and returns its old value
%
% ret = INIT( 'clear',1 )
%   Sets the values of TTEST_ERRORFLAG, TTEST_ALLTESTFLAG, TTEST_VERBOSE to their default values and clears all user defined persistent variables
%
% ret = INIT( 'v',[val] )
% ret = INIT( 'verbose',[val]
%   Sets the value of the global variable TTEST_VERBOSE to val (if val is given) and returns its old value.
%   TTEST_VERBOSE defines the verbose level of the TTEST suite
%
% INIT( 'p' )
%   Prints out the values of all global and user defined persistent variables
%
% See also: TTEST_SETUP
%
% Written by: tommsch, 2020
   
    % Install TTEST-Package
    % We autogenerate the functions, so that we do not forget any
    if( ~isequal(exist('EXPECT_EQ','file'),2) );
        TTEST_SETUP(); end;
    
    
    %Sets "environment variables" for the ttests.
    global TTEST_ERRORFLAG;         if( nargin==0 || isempty(TTEST_ERRORFLAG) );         TTEST_ERRORFLAG    = TTEST_ERRORFLAG_DEFAULT;    end;
    global TTEST_ERROR_EXPECT;      if( nargin==0 || isempty(TTEST_ERROR_EXPECT) );      TTEST_ERROR_EXPECT = TTEST_ERROR_EXPECT_DEFAULT; end;
    global TTEST_ERROR_ASSERT;      if( nargin==0 || isempty(TTEST_ERROR_ASSERT) );      TTEST_ERROR_ASSERT = TTEST_ERROR_ASSERT_DEFAULT; end;
    global TTEST_VERBOSE;           if( nargin==0 || isempty(TTEST_VERBOSE) );           TTEST_VERBOSE      = TTEST_VERBOSE_DEFAULT;      end;
    global TTEST_MAXDISP;           if( nargin==0 || isempty(TTEST_MAXDISP) );            TTEST_MAXDISP      = TTEST_MAXDISP_DEFAULT;      end;
    persistent TTEST_VAR;           if( nargin==0 || isequal(TTEST_VAR,[]) );            TTEST_VAR          = TTEST_VAR_DEFAULT;          end;
    
	
	while( ~isempty(varargin) );
        varargin{1} = lower( varargin{1} );
		switch varargin{1};
            case {'exp'};
                TTEST_ERROR_EXPECT = varargin{2};
            case {'ass'};
                TTEST_ERROR_ASSERT = varargin{2};
            case {'flag'};
                ret = TTEST_ERRORFLAG;
                if( numel(varargin)>=2 );
                    TTEST_ERRORFLAG = varargin{2}; end;
            case {'noclear'};
                %do nothing
            case {'clear'};
                TTEST_ERROR_EXPECT = TTEST_ERROR_EXPECT_DEFAULT;
                TTEST_ERROR_ASSERT = TTEST_ERROR_ASSERT_DEFAULT;
                TTEST_ERRORFLAG = TTEST_ERRORFLAG_DEFAULT;
                TTEST_VERBOSE = TTEST_VERBOSE_DEFAULT;
                TTEST_MAXDISP = TTEST_MAXDISP_DEFAULT;
                TTEST_VAR = TTEST_VAR_DEFAULT;
            case {'v','verbose'};
                ret = TTEST_VERBOSE;
                if( numel(varargin)>=2 );
                    TTEST_VERBOSE = varargin{2}; end;
            case {'p','print'};
                try;
                    fprintf( 'TTEST_ERRORFLAG      = ' ); disp( TTEST_ERRORFLAG );
                    fprintf( 'TTEST_VERBOSE        = ' ); disp( TTEST_VERBOSE );
                    fprintf( 'TTEST_ERROR_EXPECT   = ' ); disp( TTEST_ERROR_EXPECT );
                    fprintf( 'TTEST_ERROR_ASSERT   = ' ); disp( TTEST_ERROR_ASSERT );
                    fprintf( 'TTEST_MAXDISP        = ' ); disp( TTEST_MAXDISP );
                    for i = 1:numel( TTEST_VAR );
                        fprintf( '%s \t= ', TTEST_VAR{i}{1} ); disp( TTEST_VAR{i}{2} ); end;                  
                catch;
                    fprintf( 'Cannot print all global variables.\n' ); end;
            otherwise;
                found = false;
                for i = 1:numel(TTEST_VAR)
                    if( isequal(TTEST_VAR{i}{1},varargin{1}) );
                        found = true;
                        ret = TTEST_VAR{i}{2};
                        if( numel(varargin)>=2 );
                            TTEST_VAR{i}{2}=varargin{2}; end; end; end;
                if( ~found );
                    TTEST_VAR{end+1}{1} = varargin{1}; %#ok<AGROW>
                    if( numel(varargin)>=2 );
                        TTEST_VAR{end}{2} = varargin{2};
                    else
                        TTEST_VAR{end}{2} = []; end;
                    ret = []; 
                    if( TTEST_VERBOSE >= 2 );
                        fprintf( 'Create persistent variable %s\n', varargin{1} ); end; end; end;
            
        if( numel(varargin)>=2 );
            varargin(1:2) = [];
        else
            varargin(1) = []; end; end;

end

function val = TTEST_ERRORFLAG_DEFAULT;     val = false; end
function val = TTEST_ERROR_EXPECT_DEFAULT;  val = 'exp'; end
function val = TTEST_ERROR_ASSERT_DEFAULT;  val = 'ass'; end
function val = TTEST_VERBOSE_DEFAULT;       val = 1; end
function val = TTEST_MAXDISP_DEFAULT;       val = 1000; end
function val = TTEST_VAR_DEFAULT;           val = {}; end


