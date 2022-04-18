function[model] = Mallows_Ulam_learning(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = Mallows_Ulam_learning(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% Mallows_Ulam_learning:  Creates Mallows permutation model with Ulam distance
% IMPORTANT: The program assumes that the computation of the number of partitions
%            (Ferrer shapes) at each distance have been done and are saved in the files 
%            FerrerShapes_Lengths_n.mat, where n is the size of the problem.  
%            These files can be computed using the function ComputeFerrerShapes.m
%            but its computation can be costly in terms time             
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
% Created version 02/12/2015. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 04/20/2015. Roberto Santana (roberto.santana@ehu.es)  

global FitnessImprovement

initialTheta = cell2num(learning_params{1}(1)) + FitnessImprovement*cell2num(learning_params{1}(5));

upperTheta = cell2num(learning_params{1}(2));

if (initialTheta > upperTheta)
    initialTheta = upperTheta;
end


maxit = cell2num(learning_params{1}(3));

%RankingFun= cell2num(learning_params{1}(4));

NSel = size(SelPop,1); % Population size

% 0.- Check if FerreShapes are stored 
  try
    eval(['load FerrerShapes_ncounts_',num2str(NumbVar)]);   % The number of solutions at distance i will be in variable ncount
  catch
    disp(['File FerrerShapes_ncounts_',num2str(NumbVar),' is required to use this function. You can created using ComputeFerrerShapes.m ']); 
  end

% 1.- Calculate the consensus ranking/permutation, using the chosen method
%[ConsensusRanking,mindist] = eval([RankingFun,'(''ulam_distance'',SelPop,NSel)']);
[ConsensusRanking,mindist] = setMedianPermutationUlam('ulam_distance',SelPop,NSel);

% 2.- Calculate theta parameters
Theta =  CalculateThetaParameterU(ConsensusRanking,mindist,ncount,SelPop,initialTheta,upperTheta,maxit);


% 3.- Calculate Vj probabilities vector.
VProbs(1) = 1;   % exp(-theta*d) = exp(-theta *0)
for i=2:NumbVar,
   VProbs(i) = VProbs(i-1) + ncount(i)*exp(-Theta*(i-1));
end
VProbs = VProbs(2:end);
VProbs = VProbs/sum(VProbs);


model{1} = VProbs; % Vector of cumulated distance probabilities
model{2} = ConsensusRanking;
model{3} = Theta;

return;
