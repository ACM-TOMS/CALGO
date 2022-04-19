function[model] = LearnEHM(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = LearnEHM(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% LearnEHM:  Creates and Edge Histograpm Model (EHM) permutation model for symmetrical o asymmetrical
%            variants
%             
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)
% learning_params{1}(1) = Sym_VS_Asym: Specifies learning type
%                         Sym_VS_Asym=1: Symmetrical learning; Sym_VS_Asym=0: Asymmetrical learning
% learning_params{1}(2) = Bratio: Prior used for mutation-like effect (Bratio>0)

% OUTPUTS
% model: Symmetrical or Asymmetrical matrix according to type of learning
% (Sym_VS_Asym)
%
% Last version 6/23/2013. Roberto Santana (roberto.santana@ehu.es)       

Sym_VS_Asym = cell2num(learning_params{1}(1));
Bratio = cell2num(learning_params{1}(2));

N = size(SelPop,1);


EHM = zeros(NumbVar,NumbVar);
indices = [1:NumbVar,1];

% Computes the edge histogram using all solutions in the selected
% population
for i=1:N,
  for j=1:NumbVar-1,
      EHM(SelPop(i,indices(j)),SelPop(i,indices(j+1))) = EHM(SelPop(i,indices(j)),SelPop(i,indices(j+1))) + 1;
  end
end

% Computes the epsilo value according to the EHM 
epsilon = Bratio*(((Sym_VS_Asym+1)*N)/(NumbVar-1)); % Symmetrical and assymetrical EHM have different epsilons

% Applies the epsilon to the model
if Sym_VS_Asym==1    
  model{1} = EHM + EHM + epsilon*ones(NumbVar,NumbVar); 
else
  model{1} = EHM + epsilon*ones(NumbVar,NumbVar);
end


return;

