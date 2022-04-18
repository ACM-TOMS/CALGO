function [permu] = GeneratePermuFromX(x,NumbVar)
% [permu] = GeneratePermuFromX(x,NumbVar)
% GeneratePermuFromX:  Generates a random permutation from a given x vector, that contains the decomposition of the permutation.
%
% INPUTS
% x: a vector that contains,the decomposition of the Cayley distance of a permutation with the permutation identity.
%    In other words, the decomposition of the permutation.
% NumbVar:   Number of variables
% OUTPUTS
% permu: Generated permutation
%
% References:
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the cayley distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 03/27/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/03/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

inverted= -1*ones(1,NumbVar); % to know if an item is not used, when is used the value is the position in the permutation.
permu=inverted; % fill permutation with -1
randValues= rand(1,NumbVar-1); %Generate random values
%Generate a random permutation from x vector
%See Generating the permutation fordwards in [1]
for pos=1:NumbVar-1
    
    if(x(pos) == 0) % search an item that close the cycle and wasn't used
        item=pos;
        while(inverted(item) ~= -1 || ~item_closes_cycle(pos, item, permu))
            item=item-1;
        end
    else % search an item that no close the cycle and wasn't used 
        item=0;
        random=mod(randValues(pos) * intmax,NumbVar -pos);
        while(random >=0)
            item=item+1;
            if(inverted(item) == -1 && ~item_closes_cycle(pos, item, permu))
                random=random-1;
            end
        end
    end

    inverted(item)=pos;
    permu(pos)=item;
       
end
pos=pos+1;
last = find(inverted == -1); % it will be only one
permu(pos)=last;





