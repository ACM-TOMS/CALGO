function [ListFactors] =  CreateListFactorsNK(NumberVar,k)
% [ListFactors] =  CreateListFactorsNK(NumberVar,k)
% CreateListFactorsNK:   Creates the structure of an instance of the NK random landscape 
%                        The neihbors for each of the variables are
%                        randomly selected
% INPUTS
% NumbVar: Number of variables 
% k: Number of neighbors (0<k<NumbVar)
% OUTPUTS
% ListFactors: Each  cell {i} stores the current variable (i), and its k neighbors
%
% Last version 8/26/2008. Roberto Santana (roberto.santana@ehu.es)       


for i=1:NumberVar,
 PotentialNeighbors = randperm(NumberVar-1);
 %[~, PotentialNeighbors] = sort(rand(1,NumberVar-1).*rand(1,NumberVar-1));
 %e(i,:) = PotentialNeighbors;
 aux = [[1:i-1],[i+1:NumberVar]];
 PotentialNeighbors = aux(PotentialNeighbors);
 %f(i,:) = PotentialNeighbors;
 ActualNeighbors = PotentialNeighbors(1:k);
 ListFactors{i} = [i,ActualNeighbors]; 
 %d(i,:) = [i,ActualNeighbors];
 
end,

%e,f,d