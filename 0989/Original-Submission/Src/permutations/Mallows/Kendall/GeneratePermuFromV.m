function [permu] = GeneratePermuFromV(v,NumbVar)
% [permu] = GeneratePermuFromV(v,NumbVar)
% GeneratePermuFromV:  Generates the permutation corresponding to a given v vector, t
% that contains the decomposition of the permutation.
%
% INPUTS
% v: a vector that contains,the decomposition of the Kendall-tau distance of a permutation with the permutation identity.
%    In other words, the decomposition of the permutation.
% NumbVar:   Number of variables
% OUTPUTS
% permu: Generated permutation
%
% References:
% [1] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation, 2013
%
% Created version 12/21/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 03/20/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

aux_n = 1:NumbVar;

%Generate the corresponding permutation of the v vector
%Optimized version compared to the paper
for i=1:NumbVar-1

        val = v(i);
        index = 1;
        while( ~(aux_n( index ) ~= -1 && val == 0))
            if aux_n( index ) ~= -1 
                val = val-1;
			end
			index = index + 1;
		end
        permu( i ) = index  ;
        aux_n( index ) = -1 ;
end

index=1;
while(aux_n( index ) == -1 )
index = index +1;
end
permu( NumbVar) = index ;


    
