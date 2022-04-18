function [val, args, retname] = parsem( varargin )
% [val, args, retname] =  parsem( name, args , [input3], [options], ['set' || 'condset'] )
% [val, args, retname] =  parsem( args , 'test')
% Parses varargin and is easier to use than Matlabs parse.
%
% Input:
% ======
%   name                string or cell-array of strings. 
%                       The questionable argument. The function searches for string <name> in <args>. The return value depends on [input3] and on args.
%
%   args                string or cell array of strings
%         'parse_one'           If contained in args, then val = 1. Args is treated in the usual way. Thus the value in args can be different from the return value!!!
%         'parse_random'        If contained in args, then val = 0 or 1. Args is treated in the usual way. Thus the value in args can be different from the return value!!!
%         'help'                If contained in args, then 'name' is printed to the console.
%
%
%   [input3]            default value
%                           If input3 is not given: If the string can be found, and the entry in args behind the string is not a string, then this value is returned.
%                                                   Otherwise, If the string can be found the function returns 1, otherwise 0.
%                           If input3 is given:     If the string can be found, the function returns the value in <args> which is one position behind the last occurence of <name>.
%                                                   If that string cannot be found, the function returns the value which is defined input3.
%                       The last occurence of 'name' defines what happens
% Options:
% ========
%   'set'               All occurences of 'name' in args are removed and at the end 'name' (and if given input3) is appended.
%   'condset'           Only sets the value if 'name' is not found. The ordering of the element in args will be changed.
%   'test'              Tests if args is empty. If not, a warning is printed.
%                       If 'test' is given, 'name' must be ommited and no other option must be given
%   'expect',val        expected values for the <value> one position behind the last occurence of <name>
%                       If <value> is not an expected value, a warning is given
%                       only allowed if input3 is given, 
%                       format of val:
%                           {'clop',[lb ub]}             lb <= value <  ub
%                           {'opop',[lb ub]}             lb <  value <  ub
%                           {'opcl',[lb ub]}             lb <  value <= ub
%                           {'clcl',[lb ub]}             lb <= value <= ub
%                           {e1, ..., en}                value must be equal to an element in {e1,...,en} (in particular e1 must be unequal to 'clop','opop','opcl' and 'clcl')
%                           f                            f(value) must return true, where f is a function handle
%                           'name'                       value must be equal to an element in <name>. If the option <name> is not given, thus a warning is printed. 
%                                                        Thus, 'expect','name' shall only be given for mandatory options
%   'expecte',val       The same as 'expect', but an empty matrix as argument (i.e. []) is always implicitely allowed
%                       
%
% Remark:
%   If an option shall be removed, call this function with the option to be removed and use the second return argument.
%
% Output:
%   val                 The parsed value or the set value
%   args                input arguments without the parsed values (except 'set' or 'condset' is given).
%                       If the input argument is contained more than once, then all of those are removed in the return value.
%   retname             The parsed name. First entry  in 'name' if any of 'name' is not present in 'args'
%
% Eg.:  [v,arg] = parsem( 'tommsch', {'asd','tommsch',2}, 'is the best', 'set' )
%       [v,arg] = parsem( 'tommsch', {'asd','tommsch',0,'tommsch',1}, 'is the second best', 'condset' );
%       [v,arg] = parsem( 'tommsch', {'help','Olga'} );
%       parsem( 'T', {'parse_one'} );
%
% Written by: tommsch, 2016

% Changelog: tommsch, 2019-05-08    Added warning if name with parameter shall be parsed, but parameter is missing (e.g. parsem('tommsch',{23,'tommsch'},[]) )
%            tommsch, 2019-09-19    Added experimental feature: Allowing cell array of string in argument 'name'
%            tommsch, 2020-01-08    Added return value retname
%            tommsch, 2020-04-20    Behaviour change of retname
%            tommsch, 2020-04-25    Added experimental option 'expect'

if( isequal(varargin{2},'test') )
    if( ~isempty(varargin{1}) )
%         if( isequal(varargin{1},{'help'}) );
%             error( 'parsem:help', 'Early termination due to option ''help''.'); end;
        [ST,~] = dbstack;
        warning( 'parsem:unkown', 'Error in %s()/parsem(): Unkown argument(s):', ST(end).name );
        try; 
            if( isempty(varargin{1}) )
                fprintf( 'Empty argument.\n' ); 
            else
                disp(varargin{1}); end;
        catch; 
            fprintf( 'Failed to print arguments.\n' ); end;
        fprintf( '\n' ); end; %do not call vprintf, since vprintf calls parsem
    return; end;

name = varargin{1};
args = varargin{2};
if( ~iscell(args) ); 
    args = {args}; end;
if( ~iscell(name) ); 
    name = {name}; end;

idx = find( strcmp(varargin,'expect') ) ; %look, if 'allowed' is set
idxe = find( strcmp(varargin,'expecte') ) ; %look, if 'allowed' is set
if( any(idxe) );
    expecte = 1;
else
    expecte = 0; end;
idx = [idx idxe];
if( any(idx) );
    expect = varargin{idx(end)+1};
    varargin([idx idx+1]) = [];
else 
    expect = 0; end;
    

idx = find( strcmp(varargin,'set') ); %look, if 'set' is set
if( any(idx) ); 
    varargin(idx) = []; 
    set = 1; 
else; 
    set = 0; end;      

idx = find( strcmp(varargin,'condset') ); %look, if 'condset' is set
if( any(idx) ); 
    varargin(idx) = []; 
    condset = 1; 
else; 
    condset = 0; end;      

%Set default value, determine if strings are allowed for val
if( size(varargin,2)>=3 );
    input3 = varargin{3};
    %if( isa(input3,'function_handle') );
    %    fprintf( 'WARNING: Possibility to use external functions to query default values is removed.\n' ); end;
    default = input3;
    defaultgiven = 1;
else
    default = 0;
    defaultgiven = 0; end;

%read value for name, remove everything which is read 

idx = [];
if( ~isempty(name) )
    retname = name{1}; 
else
    retname = []; end;
for i = 1:numel( name );
    val = find( strcmp(args,name{i}) );
    if( ~isempty(val) );
        retname = name{i}; end;
    idx = [idx val]; end; %#ok<AGROW>
szeargs = size(args,2); %number of elements in args
if( ~isempty(idx) );
    for i = 1:size( idx, 2 );
        k = idx(i)+1;
        if( szeargs>=k && (~ischar(args{k}) || defaultgiven) );
            idx(end+1) = k;  %#ok<AGROW> %save the index of val, so that it can get deleted afterwards
            val = args{k};
        elseif( defaultgiven )
            warning( 'parsem:missing', 'In ''args'' is probably a ''value'' for some ''name'' missing.' );
            val = 1;
        else
            val = 1; end; end;
    args(idx) = [];
else
    val = default; end;

%make everything for set and condset
if( set || (condset && isempty(idx)) ); 
    args{end+1} = name{1};
    val = default;
    args{end+1} = val;
elseif( condset )
    args{end+1} = name{1};
    args{end+1} = val; end;

[parse_oneflag,parse_randomflag,helpflag] = deal( 0 );
if( any(find(strcmp(args,'parse_one'))) );
    parse_oneflag = 1;
    val = 1; end;
if( any(find(strcmp(args,'parse_random'))) );
    parse_randomflag = 1;
    val = randi( 2 )-1; end;
if( any(find(strcmp(args,'help'))) );
    helpflag = 1;
    fprintf( '       Name:   ' );
    disp( name );
    if( defaultgiven || nargout==3 )
        fprintf( '\b' );
        fprintf( '       Default: ' );
        if( nargout==3 );
           disp( name{1} );
           fprintf( '\n' );
        elseif( isempty(default) );
            if( iscell(default) );
                fprintf( '   {}\n\n' );
            else
                fprintf( '   []\n\n' ); end; 
        else
            disp( default ); end; end;
    if( ~isequal(expect,0) && ~isequal(expect,'name') );
        fprintf( '\b' );
        fprintf( '       Expect: ' );
        if( iscell(expect) && numel(expect)==2 )
            disp(expect{1}); 
            fprintf( '\b / ');
            disp(expect{2}); 
        else
            disp(expect); end; end; end;


if( ~isequal(expect,0) && ~parse_oneflag && ~parse_randomflag && ~helpflag );
    try
        if( expecte && isequal(val,[]) )
            %everything ok
        elseif( iscell(expect) && isequal(expect{1},'clop') )
            if( val<expect{2}(1) || val>= expect{2}(2) );
                warning( 'parsem:expect','Unexpected value for %s. Allowed values are %f <= val < %f', retname, expect{2}(1), expect{2}(2) ); end;
        elseif( iscell(expect) && isequal(expect{1},'opop') )
            if( val<=expect{2}(1) || val>= expect{2}(2) );
                warning( 'parsem:expect','Unexpected value for %s. Allowed values are %f < val < %f', retname, expect{2}(1), expect{2}(2) ); end;
        elseif( iscell(expect) && isequal(expect{1},'opcl') )
            if( val<=expect{2}(1) || val> expect{2}(2) );
                warning( 'parsem:expect','Unexpected value for %s. Allowed values are %f < val <= %f', retname, expect{2}(1), expect{2}(2) ); end;         
        elseif( iscell(expect) && isequal(expect{1},'clcl') )
            if( val<expect{2}(1) || val>expect{2}(2) );
                warning( 'parsem:expect','Unexpected value for %s. Allowed values are %f <= val <= %f', retname, expect{2}(1), expect{2}(2) ); end;
        elseif( iscell(expect) && ~any(cellfun(@(x)isequal(x,val),expect)) );
            s = evalc( 'disp(expect{2})' );
            warning( 'parsem:expect','Unexpected value for %s. Allowed values are%s\b', retname, s );
        elseif( isa(expect,'function_handle') && ~expect(val) )
            s = evalc( 'disp(expect)' );
            warning( 'parsem:expect','Unexpected value for %s. Allowed values are determined by%s\b', retname, s ); 
        elseif( isequal(expect,'name') && ~val );
            s = evalc( 'disp(name)' );
            warning( 'parsem:expect','Missing value from %s\b', s ); end;
    catch
        warning( 'parsem:expect','Unexpected value for %s.\n', retname ); end; end;

end

        

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   

