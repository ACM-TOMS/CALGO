function [YoungT1,YoungT2] =  GenYoungTableaux(perm)
% [YoungT1,YoungT2] =  GenYoungTableaux(perm)
% GenYoungTableaux: Given a permutation, computes the two associated Young
% Tables 
% These are auxiliary functions needed to sample new solutions of the Mallows model 
% with the Ulam distance.
%
% INPUTS
% perm: Given permutation
%
% OUTPUTS
% YoungT1: Young Table 1
% YoungT2: Young Table 2
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
% [3] J. S. Frame,  G. B. Robinson, R. M. Thrall. The hook graphs of the symmetric group. Canad. J. Math. 6, 316–325. 1954. 
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)

global YoungT1 YoungT2

%perm = [8,3,7,9,2,5,4,1,10,6];
NumbVar = length(perm);
YoungT1 = zeros(NumbVar);
YoungT2 = zeros(NumbVar);

current_col = 1;
last_row = 1;
YoungT1(1,1) = perm(1);
YoungT2(1,1) = 1;

for m=2:NumbVar,
  current_val = perm(m);  
  [current_col] =  UpdateTableaux(1,current_val,m);  
end
  
  

function [new_current_col] =  UpdateTableaux(current_col,current_val,pos)  
  global YoungT1 YoungT2
  pos_ind = find(YoungT1(:,current_col)>0);
  if isempty(pos_ind)
      YoungT1(1,current_col) = current_val;
      YoungT2(1,current_col) = pos;
      new_current_col = current_col;   
     return
  end
  sum_ind = length(pos_ind);   
  under_perm = (current_val > YoungT1(pos_ind,current_col));
  
  
  if (sum(under_perm)<sum_ind)     
    auxp = find(current_val <  YoungT1(pos_ind,current_col))
    [vmin,posmin] = min(YoungT1(auxp,current_col));    
    YoungT1(auxp(posmin),current_col) = current_val;
    %YoungT2(auxp(posmin),current_col) = pos;
    [new_current_col] =  UpdateTableaux(current_col+1,vmin,pos); 
    return
  elseif sum(under_perm)==sum_ind
     last_row = size(pos_ind,1)+1;
     YoungT1(last_row,current_col) = current_val;
     YoungT2(last_row,current_col) = pos;
     new_current_col = current_col;
  end
      
 
      
      
      
      
      
      