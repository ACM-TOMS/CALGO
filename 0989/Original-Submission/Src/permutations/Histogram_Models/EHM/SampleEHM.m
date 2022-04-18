function [NewPop] = SampleEHM(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = SampleEHM(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% SampleEHM:         Samples a population of individuals from an edge histogram model (EHM) matrix
% INPUTS
% NumbVar:   Number of variables
% model:      (model{1} = EHM matrix)
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population (May be use for partial sampling or resampling)
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% Last version 6/23/2013. Roberto Santana (roberto.santana@ehu.es)    

N = cell2num(sampling_params{1}(1)); 

% Recovers the EHM 
EHM = model{1};
NewPop=zeros(N,NumbVar);

 % The new population is generated 
 
 NewPop(:,1) = fix(rand(N,1)*NumbVar)+1; % The first variable of each solution is randomly assigned (TagNodes can be added later)
 
 
 for i=1:N,
    currentperm = NewPop(i,1); % The initial current permutation has one single city
    for j=2:NumbVar-1,         % The last city is determined by the previous (NumbVar-1)
         
          remainingperms = setdiff([1:NumbVar],currentperm); % Those cities that are not in the current permutation
          auxprob = EHM(currentperm(j-1),remainingperms); % Probabilities for remaining permutations
          auxprob = auxprob/sum(auxprob);                 % Normalization
          newpos = sus(1,cumsum(auxprob));                % Selection of the next city according to probabilies
          newpos = remainingperms(newpos);          
          currentperm = [currentperm,newpos];	          % Update of the current permutation.
    end 
    remainingperms = setdiff([1:NumbVar],currentperm);
    currentperm = [currentperm,remainingperms];           % The last city is determined by the previous ones
    NewPop(i,1:NumbVar) = currentperm;
  end            

   
 