function color = num2color(num,seed)
% color = num2color( num, seed )
% Assigns to each integer a color. For the first 7 integers, the colors are fixed, then it is pseudo-random, depending on the seed.
%
% Input:
%   num         integer
%   seed        (optional) seed for random-number generator.
%               If given, then the function does not change the state of the random-number generator.
%
% Output:
%   color       1x3 array with values in [0,1]
%
% Eg: num2color(30,2)
%
% Written by: tommsch, 2017

    switch num
        case 1; color=[1 0 0];
        case 2; color=[0 1 0];
        case 3; color=[0 0 1];
        case 4; color=[1 1 0];
        case 5; color=[1 0 1];
        case 6; color=[0 1 1];
        case 7; color=[0 0 0];
        otherwise; 
            if(nargin==1); 
                seed=1000; end;
            
            currentrng = rng;
            rng(seed); 
            rand(num);
            color=rand(1,3);
            rng(currentrng);
    end
end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.