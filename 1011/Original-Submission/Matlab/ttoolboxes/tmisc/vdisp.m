function varargout=vdisp( varargin )
% [ ret ] = vdisp( C, [str] )
% Compact (but sometimes ugly) display of 'nested arrays', 'nested structs', 'matrices', etc.
% Does not expand cell arrays inside structs.
%
% Input:
%   C                   the thing to be displayed
%   str                 (experimental) the output is appended to str and returned as ret
%
% Remarks:
%   This function is slow if the return value is used, since it calls 'evalc'.
%   This function is slow if there is a struct inside C, since it calls 'eval' and 'evalc'
%
% Eg: vdisp({[1 2 3],[1;2;3],{[1 1; 2 2; 3 3]}})
%
% See also: vprintf
%
% Written by: tommsch, 2018 
%             Stefan, University of Copenhagen



ret = ''; 
C = varargin{1};

if( nargin==2 ); 
    str = varargin{2}; 
else; 
    str = ''; end;

fs = get( 0, 'FormatSpacing' ); %store format spacing, set compact
format compact;

ret = vdisp_worker( C, ret, nargout ); %if nargout==1 we return the displayed string
if( nargout ); 
    varargout{1} = strcat( str, ret ); end;

set( 0, 'FormatSpacing',fs ); %restore format spacing

end

function ret = vdisp_worker( C, ret, saveoutput );

if( iscell(C) )
    %like celldisp, but omits linebreaks
    for i = 1:numel(C); 
        ret = vdisp_worker( C{i}, ret, saveoutput );
        if( iscell(C{1}) ); 
            if( saveoutput ); 
                %ret=[ret, sprintf('\b\n')];  %#ok<AGROW>
            else; 
                end; end; end; %fprintf('\b\n');
    return;
elseif( isstruct(C) )
    structret = dispstruct( C, [] );
    if( saveoutput ); 
        ret = [ret, structret];
    else
        disp( structret ); end;
else
    if( size(C,2)==1 && size(C,1)>1 ); 
        if( saveoutput ); 
            ret = [ret, '(transposed)']; 
        else; 
            fprintf( '(transposed)' ); end;
        C = C.'; end;
    if( saveoutput ); 
        ret = [ret, evalc( 'disp(C)' )]; 
    else; 
        disp( C ); end;
    
    if( size(C,1)>1 ); 
        if( saveoutput ); 
            ret = [ret, newline]; 
        else; 
            fprintf( '\n' ); end; end; end;
end

function fieldstring = dispstruct( structvar, substructvar )
% Display nested structures in 'structvar'; 'substructvar' is used in
%      recursion, must be [] in main call. 
% fieldstring: output
%
% Written by: Stefan, University of Copenhagen 
% Changed by: tommsch, 2019

fieldstring = [];

if isstruct( structvar )
    %display main structure at current level of recursion:
    %fieldstring=[fieldstring ' Struct ' substructvar '\n' evalc( 'disp(structvar)' )]; 
    %fieldstring=[fieldstring '\n'];
    %fieldstring=[fieldstring ' Struct ' substructvar evalc( 'disp(structvar)' )]; 
    fieldstring = [fieldstring  substructvar evalc( 'disp(structvar)' )]; 
    fieldstring = [fieldstring];
    %get fields names at current level of recursion:
    fields = fieldnames( structvar );
    for k = 1:length( fields )
        if( isstruct(eval(['structvar.' fields{k}])) )
            %recursive call to get substructures:
            %fieldstring = [fieldstring dispstruct( eval(['structvar.' fields{k}]), fields{k} )]; end; end; end;%#ok<AGROW>     
            fieldstring = [fieldstring vdisp( eval(['structvar.' fields{k}]), fields{k} )]; end; end; end;%#ok<AGROW>     
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   