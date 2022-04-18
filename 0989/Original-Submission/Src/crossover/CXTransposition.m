function [NewPop] = CXTransposition(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = CXTransposition(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% CXTransposition: 
%                              Applies a transposition operator to a population of individuals. 
%                              A transposition consists of randomly selecting a position i in the individual, an
%                              offset a, and a length l. The individual is then modified by setting 
%                              x_{i+j} = x_{i+a-j} for j={0,...,l-1}.The parameters for the transposition are received  a model 
%                              constructed with function LearnSymmetricTransposition.m
%   
% INPUTS
% NumbVar:   Number of variables
% model:     Contains four parts specifying how the transposition operator
%            is applied.  model{1} stores the indexes of individuals where
%            transposition operator will be applied. model{2} the  length of each
%            transposition. model{3} the initial location of each transposition.  
%            model{4} the offSets of each transposition
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population that it is used here to do
%            crossover between solutions
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% 
% Last version 11/7/2013. Roberto Santana (roberto.santana@ehu.es) following 
% Symmetry in evolutionary and estimation of distribution algorithms.
% (Santana, McKay, Lozano: 2013)  2013 IEEE Congress on Evolutionary
% Computation (CEC).  http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6557811&tag=1

% Input parameters are recovered
N = cell2num(sampling_params{1}(1));               % Number of solutions to sample 

TrasnIndividuals = model{1};   % Indexes of individuals where transposition operator will be applied
Len =  model{2};               % Length of each transposition
Locs = model{3};               % Location of each transposition
OffSets = model{4};            % OffSets of each transposition

DoublePop = [AuxPop,AuxPop]; % Used for garanteeing toroid structure

                                     
% For each individual the transposition operator is applied
for i=1:N,
   ind = TrasnIndividuals(i);
   loc = Locs(i);
   len = Len(i);
   offset = OffSets(i);
  % [ind,loc,len,offset]
   Individual = AuxPop(ind,:);
   for j=0:len-1,
    if (loc+j>len) 
      pos = loc+j-len;
    else
      pos = loc+j;
    end      
    Individual(pos) = DoublePop(ind,loc+offset+j);   
   end
   NewPop(i,:) = Individual;
 end  

 