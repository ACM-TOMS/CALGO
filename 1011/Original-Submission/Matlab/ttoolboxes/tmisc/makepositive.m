function M = makepositive(M); 
% M = makepositive(M); 
% Multiplies arrays such that its first non-zero entry is positive and its norm is preserved.
% i.e. if M is real, then M is multiplied by -1, if the first nonzero entry is negative.
% If M is real, then the output is also real.
%
% See also: normalizematrix
%
% E.g.: makepositive([2+1i i; 0 11])
%
% Written by: tommsch, 2018

    idx=find(M,1);
    if(isreal(M))
        if(M(idx)<0); 
            M=-M; 
        end; 
    else
        M=M*exp(-1i*angle(M(idx)));
    end
    
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   
