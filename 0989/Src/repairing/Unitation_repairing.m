function[NewPop] = Unitation_repairing(Pop,RangeValues,repairing_params)
% SetInBounds_repairing:   For a problem with binary representation
%                          Modify the solutions to the minimum (respectivel maximum) bounds 
%                          if the number of variables are under (respectively over) the variables ranges
% INPUTS
% Pop: Population 
% RangeValues: RangeValues(1) = Minimum number of Ones. RangeValues(2) = Maximum number of ones
% repairing_params: Not use for this function
% OUTPUTS
% NewPop: New repaired population
%
% Last version 8/26/2008. Roberto Santana (roberto.santana@ehu.es)       

NumbVar = size(Pop,2);
N = size(Pop,1);
NewPop = Pop;
   for i=1:N,       
          n_ones = sum(Pop(i,:));                  % Number of ones in the solution
          n_zeros = NumbVar - n_ones;              % Number of zeros in the solution
          if n_ones<RangeValues(1)                 % There are less ones than needed
            n_new_ones = RangeValues(1) - n_ones;  % How many new ones we need;  
            pos_zero = find(Pop(i,:)==0);          % Which are the current zeros                
            new_ones = randperm(n_zeros);          % Randomly ordering the zeros
            NewPop(i,pos_zero(1:n_new_ones)) = 1;      % Setting the selected zeros to one
          else n_ones>RangeValues(2)               % There are more ones than needed
            n_new_zeros = n_ones - RangeValues(2); % How many new zeros we need;  
            pos_ones = find(Pop(i,:)==1);          % Which are the current zeros                
            new_zeros = randperm(n_ones);          % Randomly ordering the zeros
            NewPop(i,pos_ones(1:n_new_zeros)) = 0;      % Setting the selected zeros to one
          end
   end   
   
% Pop = fix(2*rand(10,15))
% RangeValues = [8,8];
% [NewPop] = Unitation_repairing(Pop,RangeValues,[])
% sum(NewPop')