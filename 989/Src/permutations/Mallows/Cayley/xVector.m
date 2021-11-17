function[newvector] =  xVector(vector)
% [newvector] = xVector(vector)
% xVector:  Given a permutation xVector generates the decomposition of the permutation.
% The decomposition of the Cayley distance of a permutation with the permutation identity.
%    In other words, the decomposition of the permutation.
%
% INPUTS
% vector: a vector that contains the permutation
% OUTPUTS
% newvector: Generated decomposition of the permutation.
%
% References:
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the cayley distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 04/01/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/01/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% Calculate the decomposition of the permutation, decomposing into sum of n-1 terms, calculating the j th largest item of each cycle and assigning to newvector(j) a 0 binary value.
% See [1] for more information.

N = length(vector); %number of variables, vector length
newvector = ones(1,N);
num_cycles=0; % found number of cycles
num_visited=0;
item= 0;
visited = false(1,N); % to know item visited

while(num_visited < N)
    num_cycles = num_cycles+1;
    item=num_cycles;
    while(visited(item)) % search next non visited item from the beginning
        item=item+1;
    end

    maxItemInCycle=0;
    
    while(~visited(item)) % Calculate cycles
       if(item>maxItemInCycle)
           maxItemInCycle=item;
       end
       visited(item)=true;
       num_visited= num_visited+1;
       item = vector(item);
    end
    newvector(maxItemInCycle) = 0;
    
end

newvector=newvector(1:end-1); % the last position always is 0. The variables are defined for 1 <= j < n.

