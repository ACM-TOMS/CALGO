function x = issquare(A)
% x = issquare( A )
% Tests if an array is a (hyper)-square.
%
% E.g.: issquare(rand(3,3,3))
%
% See also: issym, iswholenumber

val=size(A);
if( all(val==val(1)) ); 
    x=true; 
else
    x=false; 
end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 