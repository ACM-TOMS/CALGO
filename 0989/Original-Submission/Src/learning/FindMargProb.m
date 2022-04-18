
function[UnivProb,BivProb]=FindMargProb(SelPop,NumbVar,Card)
% [UnivProb,BivProb]=FindMargProb(SelPop,NumbVar,Card)
% FindMargProb:      Computes the univariate and bivariate marginal probabilities
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
  UnivProb{NumbVar}(SelPop(k,NumbVar)+1) = UnivProb{NumbVar}(SelPop(k,NumbVar)+1) + 1;
  for i=1:NumbVar-1,
    UnivProb{i}(SelPop(k,i)+1) = UnivProb{i}(SelPop(k,i)+1) + 1;
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
   ijocc=iocc*jocc';  % Should be the full connection count matrix
   for i1=1:Card(i),
     for j1=1:Card(j),
       BPijInd=Card(j)*(i1-1) + j1;    %-1 + 1 to get indexing right
       BivProb{i,j}(BPijInd) = BivProb{i,j}(BPijInd)+ijocc(i1,j1);
     end;  % for j1
   end;  % for i1
  end;   % loop j
end;    % loop i
for i=1:NumbVar-1,
  for j=i+1:NumbVar,
   BivProb{i,j} = BivProb{i,j}/sum(BivProb{i,j}); % Normalization of the probabilities
 end % loop j
 UnivProb{i} = UnivProb{i}/sum(UnivProb{i});
end; % loop i
UnivProb{NumbVar} = UnivProb{NumbVar}/sum(UnivProb{NumbVar});

