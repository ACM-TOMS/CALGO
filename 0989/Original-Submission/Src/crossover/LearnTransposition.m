function [model] = LearnTransposition(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = LearnTransposition(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% LearnTransposition: Implements a GA transposition operator
%    1. randomly choosing a transpostion length len between 1 and myn/2 (i.e. half the size of the individual)
%    2. randomly choosing a location loc in the individual
%    3. randomly choosing an offset o between 2len and myn
%    4. exchanging the string starting at loc loc+len-1 with that starting at  to loc+o+len-1 (with wrapping around the toroid as necessary)
%    Learning implements Steps 1-3 (the random parameters are found each
%    individual)
%    In the sampling step the new population is created using the "model" (parameters learned from each individual) 
% INPUTS
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)
% learning_params{1} = N
%                     where:
%                     N: Number of new solutions that will be generated 
%                     BaseDimSz: Square root of the number of variables in a symmetric component
%                                used for variant 4b) of the transposition algorithm 
%                                where len is fixed to BaseDimSz^2 (i.e. the size of a basic block), but loc and o are free
%                                NOTE: When BaseDimSz=0, variant 4a (described above is applied)                                 
% OUTPUTS
% model:     Contains four parts specifying how the transposition operator
%            is applied.  model{1} stores the indexes of individuals where
%            transposition operator will be applied. model{2} the  length of each
%            transposition. model{3} the initial location of each transposition.  
%            model{4} the offSets of each transposition% Input Parameters are recovered
%
% Last version 11/7/2013. Roberto Santana (roberto.santana@ehu.es)  

N = cell2num(learning_params{1}(1));  % We need to know how many new solutions will be  generated


NSel = size(SelPop,1); % Number of selected individuals

TrasnIndividuals = fix(NSel*rand(N,1))+1;    % Indices of solutions in the Selected population where transposition will be applied

Len = fix((NumbVar/2)*rand(N,1))+1;          % Transposition length for each individual

Locs = fix(NumbVar*rand(N,1))+1;             % Location of transposition operator for each individual
DiffLenNumbVar = NumbVar - Len;              % Difference between the number of vars and length of transposition [it should be non negative]
OffSets = fix(DiffLenNumbVar .* rand(N,1))+1; % Offset for each individual
%[Len,Locs,DiffLenNumbVar,OffSets]
                                  
                                    
                                   

model{1} = TrasnIndividuals;
model{2} = Len;
model{3} = Locs;
model{4} = OffSets;






