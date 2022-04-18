function val = isAlways( val , ~, flag)
%this function exists to have a uniform interface for the function isAlways

try
    val = logical(val);
catch
    if( nargin==3 )
        switch flag
            case 'falseWithWarning'
                warning( 'double:isAlways', 'Unable to prove ''%s''.', num2str(val) );
                val = false;
            case 'false'
                val = false;
            case 'true'
                val = true;
            case 'err'
                error( 'double:isAlways', 'Unable to prove ''%s''.', num2str(val) );
            otherwise
                error( 'double:isAlways', 'Wrong Argument.' ); end;
    else   %'falseWithWarning';
        warning( 'double:isAlways', 'Unable to prove ''%s''.',num2str(val) );
        val=false; end;
    
end