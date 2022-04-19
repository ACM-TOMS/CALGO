function ret = gcdm(varargin)
% ret = gcdm( [n1,n2,n3,...,nn] )
% ret = gcdm(n1, n2, n3, ..., nn)
% Greatest common divisor of arbitrary many numbers
% Input:
%   n1,n2,...,nn        positive integers, either as a list or a vector
%
% Output
%   ret                 greatest common divisionof the input numbers
%                       gcdm([])= 0, gcdm(0,0)=0, gcdm(0,3)=3
%
% Eg: gcdm(10,15,20)
%     gcdm(0,1)
%
% See also: lcm
%
% Written by: Josh, https://de.mathworks.com/matlabcentral/profile/authors/1071611-jos

ret=0;
if(nargin==0); 
    return; end;

if(isscalar(varargin{1})); 
    ar=cell2mat(varargin); 
else
    ar=varargin{1}; 
end;

while(~isempty(ar))
    ret=gcd(ret,ar(end));
    ar(end)=[];
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 