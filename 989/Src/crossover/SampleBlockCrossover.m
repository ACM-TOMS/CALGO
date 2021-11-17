function [NewPop] = SampleSymmetricBlindBlockCrossover(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = SampleSymmetricBlindBlockCrossover(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% SampleSymmetricBlindBlockCrossover: 
%                           Applies Block-crossover where blocks are defined by equivalence classes between variables groups (e.g. symmetric components )
%                           This function assumes variables are grouped into NClasses equivalence classes of size SizClass
%                           The learning part of crossover consists of determining the indexes of 
%                           parents in the mating pools and the crossover mask for each crossover.
%                           Blind Block Crossover is similar to uniformcrossover but instead of exchanging two alleles
%                           for a gen, a complete block of variables (those that belong to the same class are exchanged between the parents) to 
%                           create the offspring
%                           A bit-flip mutation step with probability mutProb is applied after crossover

%                         TONOTE: We assume the population size is an EVEN number
% INPUTS
% NumbVar:   Number of variables
% model:     
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population that it is used here to do
%            crossover between solutions
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%

% Input parameters are recovered
N = cell2num(sampling_params{1}(1));               % Number of solutions to sample 
SymmetryIndex = cell2num(sampling_params{1}(2));   % Structure of the equivalence classes  
mutProb = cell2num(sampling_params{1}(3));         % Mutation probability [Probability of applying bit-flip to a single variable i] 

NClasses = size(SymmetryIndex,1);  % Number of classes (each class is a symmetric component)
SizClass = size(SymmetryIndex,2); % Number of elements in the symmetric component 


CrossoverMasks =  model{1}; % One masks of size NClasses for each of the N individuals
                                        % In the mask, 0 means block for first offspring taken from first parent, 1, taken from second parent 
MatingPool = model{2};  % First column, list of first parents, second column, list of second parents. They may coincide, but highly
                                     % unlikely and will produce not harm

                                     
% From each pair of parents two offspring are generated
% Blocks are exchanged according to the mask. Always within the same
% equivalence class
 for i=1:NClasses,
   indexes1 = find(CrossoverMasks(:,i)==0);
   indexes2 = find(CrossoverMasks(:,i)==1);
   NewPop(indexes1,SymmetryIndex(i,:)) = AuxPop(MatingPool(indexes1,1),SymmetryIndex(i,:));
   NewPop(N/2 + indexes1,SymmetryIndex(i,:)) = AuxPop(MatingPool(indexes1,2),SymmetryIndex(i,:));
   NewPop(indexes2,SymmetryIndex(i,:)) = AuxPop(MatingPool(indexes2,2),SymmetryIndex(i,:));
   NewPop(N/2 + indexes2,SymmetryIndex(i,:)) = AuxPop(MatingPool(indexes2,1),SymmetryIndex(i,:));
 end  

 % A bit-flip mutation step with probability mutProb is added 
 
 MutMask = find(rand(N,NumbVar) < mutProb);
 NewPop(MutMask) = 1 - NewPop(MutMask); 
