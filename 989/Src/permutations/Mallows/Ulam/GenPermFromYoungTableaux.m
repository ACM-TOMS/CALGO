function [perm] =  GenPermFromYoungTableaux(YoungT1,YoungT2)
% [perm] =  GenPermFromYoungTableaux(YoungT1,YoungT2)
% GenPermFromYoungTableaux: Given two Young Tables, generates a permutation
% according to a Mallows model that uses the Ulam distance
% A Young tableau is a Ferrers diagram with the numbers 1, 2, . . . , n in the cells such
% that each number is used once, and the entries increase along each row and down each column.
% These are auxiliary functions needed to sample new solutions of the Mallows model 
% with the Ulam distance.
%
% INPUTS
% YoungT1: Young Table 1
% YoungT2: Young Table 2
%
% OUTPUTS
% perm: Sampled permutation
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
% [3] J. S. Frame,  G. B. Robinson, R. M. Thrall. The hook graphs of the symmetric group. Canad. J. Math. 6, 316–325. 1954. 
% http://web.mit.edu/18.338/www/2012s/handouts/LIS.pdf
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)
NumbVar = length(YoungT1);
perm = zeros(1,NumbVar);

ind = 1;
while (ind<=NumbVar & sum(sum((YoungT2)))>0)
  [row,col] = find(YoungT2==max(YoungT2(:)));
  YoungT2(row,col) = 0;
  current_val = YoungT1(row,col);
  YoungT1(row,col) = 0;  
  current_col = col - 1;  
  
  while(current_col > 0)
    auxp = find(current_val(1,:) >  YoungT1(:,current_col));
    [vmax,posmax] = max(YoungT1(auxp,current_col));           
    YoungT1(posmax,current_col) = current_val;
    current_val = vmax;
    current_col = current_col - 1;
  end
  perm(ind) = current_val;    
  ind = ind + 1;      
end

perm = perm(end:-1:1);

      
      
      
      
      
      