function x = iswholenumber(A)
% [ x ] = iswholenumber( A )
% Tests if an array contains only whole numbers.
% x=~isinf(A) & floor(A)==A;
%
% Written by: tommsch, 2018

if(issym(A))
    try;
        A=double(A);
    catch;
        x=false;
        return;
    end;
end
if(~isnumeric(A)); 
    x=false; 
    return; 
else;
    x=~isinf(A) & floor(A)==A;
end;

end
