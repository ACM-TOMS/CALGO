function[model] = LearnNHM(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = LearnNHM(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% LearnNHM:  Creates NHM permutation model for symmetrical o asymmetrical
%            variants
%             
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)

% learning_params{1}(1) = Bratio: Prior used for mutation-like effect (Bratio>0)

% OUTPUTS
% model: Symmetrical or Asymmetrical matrix according to type of learning
% (Sym_VS_Asym)
%
% Last version 8/26/2008. Roberto Santana (roberto.santana@ehu.es)       


Bratio = cell2num(learning_params{1}(1));
N = size(SelPop,1);

% Computation of the histogram model;
NHM = hist(SelPop,NumbVar)';

epsilon = Bratio*(N/NumbVar); 
model{1} = NHM + epsilon*ones(NumbVar,NumbVar);


return;

