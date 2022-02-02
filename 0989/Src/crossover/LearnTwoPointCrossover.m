function [model] = LearnTwoPointCrossover(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = LearnTwoPointCrossover(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% LearnTwoPointCrossover: Determines the MatingPool and crossover points of
%                         a crossover algorithm. 
% INPUTS
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)
% learning_params{1} = SymmetryIndex;
% OUTPUTS
% model:      model{1} contains a matrix with one row for each individual, and two columns.
%             The first column is the list of  the first parents. The second column, the list of second parents. 
%             They may coincide, but this is highly unlikely and will produce not harm
%             model{2} contains,  for each pair of solutions, the point
%             where crossover is applied to
%

% Input Parameters are recovered
N = cell2num(learning_params{1}(1));  % We need to know how many new solutions will be  generated


NSel = size(SelPop,1); % Number of selected individuals

MatingPool = fix(NSel*rand(N/2,2))+1;  % First column, list of first parents, second column, list of second parents. They may coincide, but highly
                                     % unlikely and will produce not harm
Points(:,1) = fix(NumbVar*rand(N/2,1))+1;                                     
Points(:,2) = Points(:,1)+fix((NumbVar-Points(:,1)).*rand(N/2,1));
                                   

model{1} = MatingPool;
model{2} = Points;






