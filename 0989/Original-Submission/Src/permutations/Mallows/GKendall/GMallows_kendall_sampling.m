function [NewPop] = GMallows_kendall_sampling(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = GMallows_kendall_sampling(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% GMallows_kendall_sampling:  Samples a population of individuals from an Generalized Mallows kendal model
%             
% INPUTS
% NumbVar:   Number of variables
% model:      model{1} contains a permutation probabilities matrix with NumbVar-1 rows and NumbVar columns.
%					   The matrix is a left upper triangular matrix.
%             model{2} contains the consensus ranking.
%             model{3} contains calculated theta value.
%             model{4} contains a vector with calculated psi values.
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population (May be use for partial sampling or resampling)
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% References:
% [1] C. L. Mallows: Non-null ranking models. Biometrika, 1957
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation, 2013
%
% Created version 04/11/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 


VProbs = model{1}; % Permutation probabilities matrix
ConsensusRanking = model{2};
Theta = model{3};
Psis = model{4};

N = cell2num(sampling_params{1}(1)); 

%sampling stage auxiliary variables
aux = -1*ones(1,NumbVar);
aux_v = zeros(1,NumbVar); %last position always 0.

NewPop=zeros(N,NumbVar);

randValues= rand(N,NumbVar-1); %Generate random values
limit = [NumbVar-1:-1:1]; %Calculate limit vector for later use
for j=1:N,
	% generate samples and calcualte likelihood

	for i=1:NumbVar-1
		% sampling of n-1 positions of the V vector that define the permutation
		randVal=randValues(j,i); %Obtain random value
		acumul=VProbs(i,1);
		
		index=0;
		
		while (index < limit(i)) && (acumul < randVal)
			acumul = acumul + VProbs(i,(index+1)+1);
			index = index + 1;
		end
		aux_v(i)=index;
		
	end

	aux_v(NumbVar) = 0; %last position always is 0.

	permu = GeneratePermuFromV(aux_v,NumbVar); % Generate new permutation from V vector

	% Calculate composition
	newpermutation = Compose(permu,ConsensusRanking);
	
	NewPop(j,1:NumbVar) = newpermutation;
end

return;
