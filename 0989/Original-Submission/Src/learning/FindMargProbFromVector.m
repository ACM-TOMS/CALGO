
function[UnivProb,BivProb]=FindMargProbFromVector(SelPop,NumbVar,Card,Vector)
% [UnivProb,BivProb]=FindMargProbFromVector(SelPop,NumbVar,Card,Vector)
% FindMargProb:      Computes the univariate and bivariate marginal
%                    probabilities from a vector of probabilities (e.g. Boltzmann)
%                    associated to the selected population 
% INPUT:
% SelPop: Selected population 
% NumbVar: Number of variables
% Card: Vector with the dimension of all the variables. 
% OUTPUT:
% UnivProb:  Univariate probabilities
% BivProb:  Bivariate probabilities
% Last version 12/13/2008. Bob Mckay   


SelPopSize = size(SelPop,1);
for i=1:NumbVar-1,
 UnivProb{i} = zeros(1,Card(i));   %Initialization of univariate probabilities
 if (i==NumbVar-1)
    UnivProb{NumbVar} = zeros(1,Card(NumbVar));
 end; % if i
 for j=i+1:NumbVar,
   BivProb{i,j} = zeros(1,Card(i)*Card(j));  %Initialization of bivariate probabilities for all pair of variables are computed
  end; %loop j
end; %loop i
for k=1:SelPopSize,
  UnivProb{NumbVar}(SelPop(k,NumbVar)+1) = UnivProb{NumbVar}(SelPop(k,NumbVar)+1) + Vector(k);
  for i=1:NumbVar-1,
    UnivProb{i}(SelPop(k,i)+1) = UnivProb{i}(SelPop(k,i)+1) + Vector(k);
    SelPopki=SelPop(k,i);
  end;   % for i
end;     % for k

for i=1:NumbVar-1,
  iocc=zeros(Card(i),SelPopSize);
  for j=i+1:NumbVar,
   jocc=zeros(Card(j),SelPopSize);
   ijocc=zeros(Card(i),Card(j));
   for i1=1:Card(i),
     iocc(i1,1:SelPopSize)=(SelPop(:,i)' == (i1-1)*ones(1,SelPopSize));
   end;  % for i1
   for j1=1:Card(j),
     jocc(j1,1:SelPopSize)=(SelPop(:,j)' == (j1-1)*ones(1,SelPopSize));
   end;  % for j1
   
   for i1=1:Card(i),
     for j1=1:Card(j),       
       CommonPos = find((iocc(i1,:)+jocc(j1,:))==2);  % Positions that satisfy the marginal
       ijocc(i1,j1) = sum(Vector(CommonPos));     
       BPijInd=Card(j)*(i1-1) + j1;    %-1 + 1 to get indexing right
       BivProb{i,j}(BPijInd) = BivProb{i,j}(BPijInd)+ijocc(i1,j1);
     end;  % for j1
   end;  % for i1
  end;   % loop j
end;    % loop i
alpha = 0.8;
for i=1:NumbVar-1,
  for j=i+1:NumbVar,      
    BivProb{i,j} = alpha*BivProb{i,j} + (1-alpha)*(1/9);   
    BivProb{i,j} = BivProb{i,j}/sum(BivProb{i,j}); % Normalization of the probabilities
    %  BivProb{i,j} 
  end % loop j
  UnivProb{i} = alpha*UnivProb{i} + (1-alpha)*(1/3);
  UnivProb{i} = UnivProb{i}/sum(UnivProb{i});
  %UnivProb{i}
end; % loop i
UnivProb{NumbVar} = alpha*UnivProb{NumbVar} + (1-alpha)*1/3;
UnivProb{NumbVar} = UnivProb{NumbVar}/sum(UnivProb{NumbVar});

