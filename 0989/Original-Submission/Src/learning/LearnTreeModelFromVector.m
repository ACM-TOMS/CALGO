function [model] = LearnTreeModelFromVector(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% [model] = LearnTreeModelFromVector(k,NumbVar,Card,SelPop,AuxFunVal,learning_params)
% LearnTreeModelFromVector: A maximum weighted tree is learned from the matrix of mutual information between the variables
%                           from a vector of probabilities (e.g. Boltzmann) associated to the selected population 
%
% INPUTS
% k: Current generation
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% SelPop:  Population from which the model is learned 
% AuxFunVal: Evaluation of the data set (required for some learning algorithms, not for this one)
% learning_params{1}(1) = {Vector}
% OUTPUTS
% model: Structure containing the tree structure (model{1} = Cliques)
%        and the parameters (model{2} = Tables)
%
% Last version 5/6/2009. Roberto Santana (roberto.santana@ehu.es)       

global BoltzMannProb;

Vector = BoltzMannProb;
  N = size(SelPop,1);
  Cliques =  [];   

      % Univariate and Bivariate probabilities are learned   
        [UnivProb,BivProb]= FindMargProbFromVector(SelPop,NumbVar,Card,Vector);
       
       % The Matrix of Mutual information is learned
         MI = IntMutualInfFromMargProb(NumbVar,Card,UnivProb,BivProb);
         
       % The tree is extracted from the matrix of mutual information
         [Cliques] = CreateTreeStructure(MI);       
         
      for j=1:NumbVar        
        nparent = Cliques(j,1);        
        if  nparent== 0 % The variable has no parent
         i = Cliques(j,3);
         Tables{j} = (UnivProb{i}*N + 1) ./ ( (N + Card(i))*ones(1,Card(i)));
         %UnivProb{i} = Tables{i};
        else 
         i = Cliques(j,4);
         parent = Cliques(j,3);
         %[i,parent,Card(i),Card(parent)]            
         if parent<i
           a = BivProb{parent,i};
           %AuxBivProb = reshape(a',Card(parent),Card(i))
           AuxBivProb = reshape(a',Card(i),Card(parent))';
         else
           a = BivProb{i,parent};
           AuxBivProb = reshape(a',Card(parent),Card(i));           
         end
            aux = repmat(UnivProb{parent}',1,Card(i));
            %auxcard = repmat(Card(i),Card(parent),Card(i));     
            %CondBivProb =  (AuxBivProb*N +1) ./ (aux*N + auxcard);  % Laplace Estimator     
            LapBivProb  = (AuxBivProb*N +1)/(N+Card(i)*Card(parent));           
            CondBivProb =  LapBivProb ./ aux;  % Laplace Estimator        
            CondBivProb =  CondBivProb ./ repmat(sum(CondBivProb')',1,Card(i));
            CondBivProb = AuxBivProb ./aux;
            Tables{j} = CondBivProb;            
        end        
      end
      
     
     model{1} = Cliques;
     model{2} = Tables;
     
     return;
     
     
         