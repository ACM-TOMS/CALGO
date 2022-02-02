function[SelPop,SelFunVal]= Boltzmann_selection(Pop,FunVal,selection_params)
% [SelPop,SelFunVal]= prop_selection(Pop,FunVal,selection_params)
%
% prop_selection:       Proportional selection for single objective functions
% INPUTS 
% Pop:                          Original population
% FunVal:                       A matrix of function evaluations, one vector of m objectives for each individual
% OUTPUTS
% SelPop: Selected population
% SelFunVal:  A vector of function evaluations for each selected individual
%
% Last version 8/26/2008. Roberto Santana (roberto.santana@ehu.es)       
 
global BoltzMannProb

   T = cell2num(selection_params{1}(1));
   PopSize = size(Pop,1);
   
   SelPop = Pop;
   SelFunVal = FunVal;
   
  
   if min(FunVal)==max(FunVal)
     BoltzMannProb = 1/PopSize*ones(PopSize,1);  
   else
     NormFunVal = (FunVal-min(FunVal))/(max(FunVal)-min(FunVal));
     
     BoltzMannProb = exp(NormFunVal/T);   
     BoltzMannProb = BoltzMannProb/sum(BoltzMannProb);
   end
  
      
   return
 
 
 