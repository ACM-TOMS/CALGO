function [NewPop] = Mallows_cayley_sampling(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = Mallows_cayley_sampling(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% Mallows_cayley_sampling:  Samples a population of individuals from an Mallowa cayley model
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
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 


XProbs = model{1}; % Permutation probabilities matrix
ConsensusRanking = model{2};
Theta = model{3};
Psis = model{4};

N = cell2num(sampling_params{1}(1)); 

%sampling stage auxiliary variables
aux = -1*ones(1,NumbVar);
aux_x = zeros(1,NumbVar); %last position always 0.

NewPop=zeros(N,NumbVar);

randValues= rand(N,NumbVar-1); %Generate random values

for j=1:N,
	% generate samples and calcualte likelihood

    i=1:NumbVar-1;
    % sampling of n-1 positions of the V vector that define the permutation
    randVal=randValues(j,i); %Obtain random values
    
    aux_x=(randVal >= XProbs(i)'); %if randVal < XProbs(i) then aux_x(i) = 0;

	aux_x(NumbVar) = 0; %last position always is 0.

	permu = GeneratePermuFromX(aux_x,NumbVar); % Generate new permutation from X vector

	% Calculate composition
	newpermutation = Compose(permu,ConsensusRanking);
	
	NewPop(j,1:NumbVar) = newpermutation;
end

return;
