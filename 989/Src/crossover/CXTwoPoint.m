function [NewPop] = CXTwoPoint(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = CXTwoPoint(NumbVar,model,Card,AuxPop,AuxFunVal,CX_params)
% CXTwoPoint: Implements the crossover part of LearnTwoPointCrossover
%             
%               
% INPUTS
% NumbVar:   Number of variables
% model:      model{1} contains a matrix with one row for each individual, and two columns.
%             The first column is the list of  the first parents. The second column, the list of second parents. 
%             They may coincide, but this is highly unlikely and will produce not harm
%             model{2} contains,  for each pair of solutions, the point where crossover is applied to
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


MatingPool = model{1};  % First column, list of first parents, second column, list of second parents. They may coincide, but highly
                                     % unlikely and will produce not harm
Points = model{2};      % For each pair of solutions, the point where crossover is applied to

% From each pair of parents two offspring are generated
for i=1:N/2,
   % First offspring is generated 
   %[i,Points(i,:)]
   NewPop(i,1:Points(i,1)) = AuxPop(MatingPool(i,1),1:Points(i,1)); % First segment taken from Parent 1
   NewPop(i,Points(i,1)+1:Points(i,2)) = AuxPop(MatingPool(i,2),Points(i,1)+1:Points(i,2)); % Second segment taken from Parent 2   
   NewPop(i,Points(i,2)+1:NumbVar) = AuxPop(MatingPool(i,1),Points(i,2)+1:NumbVar); % Third segment taken from Parent 1
   
   
   % Second offspring is generated
   NewPop(i+N/2,1:Points(i,1)) = AuxPop(MatingPool(i,2),1:Points(i,1)); % First segment taken from Parent 2
   NewPop(i+N/2,Points(i,1)+1:Points(i,2)) = AuxPop(MatingPool(i,1),Points(i,1)+1:Points(i,2)); % Second segment taken from Parent 1
   NewPop(i+N/2,Points(i,2)+1:NumbVar) = AuxPop(MatingPool(i,2),Points(i,2)+1:NumbVar); % Third segment taken from Parent 2
end  

 
