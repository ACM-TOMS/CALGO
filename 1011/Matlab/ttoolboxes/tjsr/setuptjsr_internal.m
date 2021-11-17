function setuptjsr_internal(flag)
% setuptjsr_internal( [flag] )
% Tests some functions used by tjsr()
% 
% Options:
%       flag        if set, makes a full test
%                   A higher value for flag, also increases the verbose level
%
% Written by: tommsch, 2019

%#ok<*NOPRT>
%#ok<*ASGLU>
%#ok<*NASGU>
 %#ok<*ALIGN>
 
if(nargin==0); flag=0; end;
selftesttjsr_internal(flag);

end

function selftesttjsr_internal(flag)

    fprintf('Some tjsr_ - functions will be tested: ');
    if(~isequal(TJSR_CONEFUNCT,0)); vprintf('Wrong value for TJSR_CONEFUNCT.\n'); end;
    if(~isequal(TJSR_MINKFUNCT,1)); vprintf('Wrong value for TJSR_MINKFUNCT.\n'); end;
    if(~isequal(TJSR_COMPLEXFUNCT,2)); vprintf('Wrong value for TJSR_COMPLEXFUNCT.\n'); end;
    if(~isequal(tjsr_zeroJsr({[0 1 0;0 0 1;0 0 0],[0 0 2;0 0 3;0 0 0]}),true));  vprintf('tjsr_zeroJsr broken.\n'); end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 