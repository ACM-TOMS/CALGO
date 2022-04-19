function[SelPop,SelFunVal]=NonDominated_selection(Pop,FunVal,selection_params)
% [SelPop,SelFunVal]=NonDominated_selection(Pop,FunVal,selection_params)
% NonDominated_selection:      The set of non dominated individuals is selected.                          
%                              (in some cases this set could be very small
% INPUTS 
% Pop:                 Original population
% FunVal:              A matrix of function evaluations, one vector of m objectives for each individual
% selection_params{1}: Truncation parameter T \in (0,1)
% OUTPUTS
% SelPop: Selected population
% SelFunVal:  A vector of function evaluations for each selected individual
%
% Last version 8/26/2008. Roberto Santana (roberto.santana@ehu.es)       
 
   
   [Index]=FindParetoSet(Pop,FunVal)
   
   SelPop = Pop(Index,:);  
   SelFunVal = FunVal(Index,:); 
  
   return
 
   
   
   