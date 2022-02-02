function[model] = Mallows_kendall_learning(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = Mallows_kendall_learning(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% Mallows_kendall_learning:  Creates Mallows permutation model with kendall distance
%             
% INPUTS
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)
% learning_params{1}(1) = initialTheta: Initial theta value
% learning_params{1}(2) = upperTheta: The maximum theta value
% learning_params{1}(3) = maxit: Maximum iteration
% learning_params{1}(4) = RankingFun: The ranking function
% learning_params{1}(5) = DeltaInitialTheta

% OUTPUTS
% model:      model{1} contains a permutation probabilities matrix with NumbVar-1 rows and NumbVar columns.
%					   The matrix is a left upper triangular matrix.
%             model{2} contains the consensus ranking.
%             model{3} contains calculated theta value.
%             model{4} contains a vector with calculated psi values.
%
% References:
% [1] C. L. Mallows: Non-null ranking models. Biometrika, 1957
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
%
% Created version 12/19/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

global FitnessImprovement

initialTheta = cell2num(learning_params{1}(1)) + FitnessImprovement*cell2num(learning_params{1}(5));

upperTheta = cell2num(learning_params{1}(2));

if (initialTheta > upperTheta)
    initialTheta = upperTheta;
end

maxit = cell2num(learning_params{1}(3));

RankingFun= char(cellstr(learning_params{1}(4)));

NSel = size(SelPop,1); % Population size

% 1.- Calculate the consensus ranking/permutation, using the chose method
ConsensusRanking = eval([RankingFun,'(''kendall_distance'',SelPop,NSel)']);

% 2.- Calculate theta parameters
Theta =  CalculateThetaParameterK(ConsensusRanking,SelPop,initialTheta,upperTheta,maxit);

% 3.- Calculate psi constants.
Psis = CalculatePsiConstantsK(Theta, NSel);

% 4.- Calculate Vj probabilities matrix.
VProbs = VProbMat(NumbVar,Psis,Theta);

model{1} = VProbs; % Permutation probabilities matrix
model{2} = ConsensusRanking;
model{3} = Theta;
model{4} = Psis;


return;
