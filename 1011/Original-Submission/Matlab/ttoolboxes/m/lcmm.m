function ret = lcmm(varargin)
% ret = lcmm( [n1,n2,n3,...,nn] )
% ret = lcmm(n1, n2, n3, ..., nn)
% Least common multiple of arbitrary many numbers
% Input:
%   n1,n2,...,nn        positive integers, either as a list or a vector
%
% Output
%   ret                 least common multiple of the input numbers
%                       lcmm([])= 1, lcmm(0,10)=1, lcmm(-4)=4,
%
% Eg: lcmm(10,15,20)
%     lcmm(0,1)
%
% See also: lcm
%
% Written by: Josh, https://de.mathworks.com/matlabcentral/profile/authors/1071611-jos

 %#ok<*ALIGN>

    if(nargin==0); 
        ret=1; 
        return; end;

    if(isscalar(varargin{1})); 
        numberArray=cell2mat(varargin); 
    else
        numberArray=varargin{1}; 
    end;
    numberArray=abs(numberArray);

    numberArray = reshape(numberArray, 1, []);

    if(size(numberArray,2)==0); 
        ret=1; 
        return; end;

    if(any(numberArray==0)); 
        ret=0;
        return; end;

    % prime factorization array
    for i = 1:size(numberArray,2)
        temp = factor(numberArray(i));

        for j = 1:size(temp,2)
            ret(i,j) = temp(1,j); %#ok<AGROW>
        end
    end

    % generate prime number list
    p = primes(max(max(ret)));
    % prepare list of occurences of each prime number
    q = zeros(size(p));

    % generate the list of the maximum occurences of each prime number
    for i = 1:size(p,2)
        for j = 1:size(ret,1)
            temp = length(find(ret(j,:) == p(i)));
            if(temp > q(1,i))
                q(1,i) = temp;
            end
        end
    end

    % the algorithm
    z = p.^q;

    ret = 1;

    for i = 1:size(z,2)
        ret = ret*z(1,i);
    end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 