function[model] = GMallows_kendall_learning(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = GMallows_kendall_learning(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% GMallows_kendall_learning:  Creates Generalized Mallows permutation model with kendall distance
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
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation, 2013
%
% Created version 04/11/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

global FitnessImprovement

initialTheta = cell2num(learning_params{1}(1)) + FitnessImprovement*cell2num(learning_params{1}(5));

upperTheta = cell2num(learning_params{1}(2));

if (initialTheta > upperTheta)
    initialTheta = upperTheta;
end

maxit = cell2num(learning_params{1}(3));

RankingFun= char(cellstr(learning_params{1}(4)));

NSel = size(SelPop,1); % Population size

% 1.- Calculate the consensus ranking/permutation, using the chosen method
ConsensusRanking = eval([RankingFun,'(''kendall_distance'',SelPop,NSel)']);

% 2.- Calculate theta parameters
Theta =  CalculateThetaParameterGK(ConsensusRanking,SelPop,initialTheta,upperTheta,maxit);

% 3.- Calculate psi constants.
Psis = CalculatePsiConstantsGK(Theta, NSel);

% 4.- Calculate Vj probabilities matrix.
VProbs = VProbMat(NumbVar,Psis,Theta);

model{1} = VProbs; % Permutation probabilities matrix
model{2} = ConsensusRanking;
model{3} = Theta;
model{4} = Psis;


return;
