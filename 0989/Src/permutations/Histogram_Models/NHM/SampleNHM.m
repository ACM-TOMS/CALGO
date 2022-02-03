function [NewPop] = SampleNHM(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = SampleNHM(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% SampleNHM:         Samples a population of individuals from an NHM matrix
% INPUTS
% NumbVar:   Number of variables
% model:      (model{1} = NHM matrix)
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population (May be use for partial sampling or resampling)
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% Last version 8/26/2008. Roberto Santana (roberto.santana@ehu.es)    

N = cell2num(sampling_params{1}(1)); 

NHM = model{1};
NewPop=zeros(N,NumbVar);

 % The new population is generated
 
 for i=1:N,
    currentperm = [];           % The initial current permutation is empty
    for j=1:NumbVar-1,         %  Only NumbVar-1 cities, the last city is determined by the previous
         
          remainingperms = setdiff([1:NumbVar],currentperm); % Those cities that are not in the current permutation
          auxprob = NHM(j,remainingperms); % Probabilities for remaining permutations
          auxprob = auxprob/sum(auxprob);                 % Normalization
          newpos = sus(1,cumsum(auxprob));                % Selection of the next according to probabilies
          newpos = remainingperms(newpos);          
          currentperm = [currentperm,newpos];	          % Update of the current permutation.
    end 
    remainingperms = setdiff([1:NumbVar],currentperm);
    currentperm = [currentperm,remainingperms];           % The last city is determined by the previous
    NewPop(i,1:NumbVar) = currentperm;
  end            

   
 