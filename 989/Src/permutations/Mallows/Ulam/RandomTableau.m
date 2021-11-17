function [FinalTableau] = RandomTableau(NumbVar,shape)
% [Tableau] = RandomTableau(NumbVar,shape)
% RandomTableau(NumbVar,shape):  Samples a random Young Table among those of a given Ferrer diagram shape.
% 
% INPUTS
% NumbVar: Number of variables
% shape: Shape of the Tableau
% OUTPUTS
% Tableau: Sampled Tableau
% 
% References:
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the cayley distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
% [3] J. S. Frame,  G. B. Robinson, R. M. Thrall. The hook graphs of the   symmetric group. Canad. J. Math. 6, 316–325. 1954. 
%
% Created version 31/05/2016. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 31/05/2016. Roberto Santana (roberto.santana@ehu.es)  

 Tableau(1:NumbVar) = 0;       % Tableau will contains na unique numbers where na is the number of
                               % rows in the Ferrer shape
  count=0;                     % If Tableau(i) = a it means that number i is included in row a 
  for i=1:length(shape),    
    for j=1:shape(i),    
      Tableau(j) = Tableau(j) + 1;
      count = count + 1;
    end
  end   

  
  for l=1:NumbVar,
    notfound = 0;  
    while(notfound==0)
      i = randi(Tableau(1));
      j = randi(shape(1));      
      notfound = (i<= Tableau(j) & j <= shape(i));       
    end
    
    
    ih = 1;
    while(ih~=0)
      ih = Tableau(j) + shape(i) - i - j;
      if(ih>0)
        k = randi(ih);
        if ( k <= shape(i)-j)
          j = j + k;
        else
          i = k - shape(i) + i + j;
        end
      end  
    end
    shape(i) = shape(i) - 1;
    Tableau(j) = Tableau(j) - 1;
    Tableau(NumbVar+1-l) = i;
  end
  
  nvals = unique(Tableau);
  nrows = length(nvals);
  
  FinalTableau = zeros(nrows,NumbVar);    % FinalTableau converts to the Ferrer shape structure
  for i=1:nrows,                          % It is only a different representation of Tableau
    aux = find(Tableau==nvals(i));
    FinalTableau(i,1:length(aux)) = aux;
  end   

  return
