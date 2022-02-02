function [model] = LearnSymmetricBlindBlockCrossover(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = LearnSymmetricBlindBlockCrossover(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% LearnEnsembleTreeModel: This function assumes variables are grouped into  NClasses equivalence classes of size SizClass
%                         The learning part of crossover consists of determining the indexes of
%                         parents in the mating pools and the crossover mask for each crossover.
%                         Blind Block Crossover is similar to uniformcrossover but instead of exchanging two alleles
%                         for a gen, a complete block of variables (those that belong to the same class are exchanged between the parents) to 
%                         create the offspring
%                         TONOTE: We assume the population size is an EVEN number
% INPUTS
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)
% learning_params{1} = SymmetryIndex;
% OUTPUTS
% model: Set of (NClasses) tree structures model = EnsembleOfModels
%        where  (EnsembleOfModels{i}{1} = Cliques) contains the tree
%        structure  and (EnsembleOfModels{i}{2} = Tables)
%        the tree parameters.
%

% Input Parameters are recovered
SymmetryIndex = cell2num(learning_params{1}(1));
N = cell2num(learning_params{1}(2));  % We need to know how many new solutions will be  generated

NClasses = size(SymmetryIndex,1);  % Number of classes (each class is a symmetric component)
SizClass = size(SymmetryIndex,2);  % Number of elements in the symmetric component 


NSel = size(SelPop,1); % Number of selected individuals

CrossoverMasks = fix(rand(N/2,NClasses)); % One masks of size NClasses for each of the N individuals
                                        % In the mask, 0 means block for first offspring taken from first parent, 1, taken from second parent 
MatingPool = fix(NSel*rand(N/2,2))+1;  % First column, list of first parents, second column, list of second parents. They may coincide, but highly
                                     % unlikely and will produce not harm
                                     
                                    
                                   

model{1} = CrossoverMasks;
model{2} = MatingPool;






