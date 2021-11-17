function[Tables] = LearnFDAParameters(Cliques,SelPop,NumbVar,Card)
% [Tables] = LearnFDAParameters(Cliques,SelPop,NumbVar,Card)
% LearnFDAParameter: The parameters of the factorized  model represented in Cliques 
%                    are learned from data
% INPUTS
% Cliques: Structure of the model in a list of cliques that defines the junction graph. 
%          Each row of Cliques is a clique. The first value is the number of overlapping variables. 
%          The second, is the number of new variables.
%          Then, overlapping variables are listed and  finally new variables are listed.
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned         
% OUTPUTS
% Tables: Probability tables for each variable conditioned on its neighbors
%
% Last version 8/26/2008. Roberto Santana and Siddarta Shakya (roberto.santana@ehu.es)    


NumberCliques = size(Cliques,1);


%%%%%%%%%%%%%%%%%%%%%%  First step  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The tables of all the cliques are filled  


for i=1:size(Cliques,1),

  sizeCliqOther = Cliques(i,2);
  sizeCliqSolap = Cliques(i,1);

  CliqOther = Cliques(i,sizeCliqSolap+3:sizeCliqSolap+sizeCliqOther+2);
  AccCardOther = FindAccCard(sizeCliqOther,Card(CliqOther));
  dimOther =   NumconvertCard(Card(CliqOther)-1,sizeCliqOther,AccCardOther)+1;

  if(sizeCliqSolap > 0)
    CliqSolap = Cliques(i,3:(sizeCliqSolap+2));
    AccCardSolap = FindAccCard(sizeCliqSolap,Card(CliqSolap));
    dimSolap =   NumconvertCard(Card(CliqSolap)-1,sizeCliqSolap,AccCardSolap)+1;
    aux=zeros(dimSolap,dimOther);
  else 
    AccCardSolap = [];
    CliqSolap = [];
    aux=zeros(1,dimOther);
    dimSolap = 1;
  end
   
  AllVars = [CliqSolap,CliqOther];

 
 for j=1:dimSolap
    if (sizeCliqSolap>0) 
      solapval = IndexconvertCard(j-1,sizeCliqSolap,AccCardSolap);
    else
      solapval=[];
    end

    for k=1:dimOther
     auxSelPop=SelPop(:,[CliqSolap,CliqOther]);
     otherval =  IndexconvertCard(k-1,sizeCliqOther,AccCardOther);
     allvarvalues = [solapval,otherval];
    
     if(size(allvarvalues,2)==1)
       aux(j,k) =  sum((auxSelPop==repmat(allvarvalues,size(SelPop,1),1))');
     else 
       aux(j,k)=sum( sum((auxSelPop==repmat(allvarvalues,size(SelPop,1),1))') == size(allvarvalues,2)); 
     end
    end 
   %aux(j,:) = (aux(j,:))/(sum(aux(j,:)));
   aux(j,:) = (aux(j,:)+1)/(sum(aux(j,:))+dimOther); % Laplace Estimator is applied together with the computation of marginal probabilities
  end
 %aux=aux/sum(sum(aux)); % Normalization
 
 
 % In Table i the probabilities of clique i are stored

 eval(['Tables{',num2str(i),'}=aux;']);
         
end

