function [objPoly,ineqPolySys,lbd,ubd] ...
    = randomUnconst(nDim,minCliqueSize,maxCliqueSize,maxDegree,randSeed); 

% [objPoly,ineqPolySys,lbd,ubd]=
% randomUnconst(nDim,minCliqueSize,maxCliqueSize,maxDegree,randSeed);
% 
% A randomly generated unconstrained problem, which is
% described in "Sums of Squares and Semidefinite Programming
% Relaxations for Polynomial Optimization Problems with Structured
% Sparsity", B-411
%
% <Input> 
% nDim: This is the dimension of the function
% minCliqueSize: minimum value of the size of the cliques generated
%                by function "chodalStruct" in unconProb.m.
%
% maxCliqueSize: maximum value of the size of the cliques generated
%                by function "chodalStruct" in unconProb.m.
%
% maxDegree: the degree of objective function, this number is even.
%
% randSeed: the seed to generate the random numbers.
%
% <Output>
% objPoly,ineqPolySys,lbd,ubd
%
% This function outputs objective functin of the following
% unconstrained polynomial programming:
% 
% min \sum_{i=1}^{noOfCliques}\{\v_i(\x)^T\V_i\v_i(\x) + g_i(\x)+\}
% sub.to \x in \Real^{nDim},
%
% where \v_i(\x) contains x_j^{maxDegree/2} (j \in Clique_i), 
% g_i(\x) is sparse polynomial and has only monomials that
% constructed by x_j (j \in Clique_i), deg g_i < maxDegree, 
% and \V_i is psd matrix.
%
% This function is used in 'solveRandUnconst.m'.
%
%
% <Input value>
% nDim:           the number of variable of objective function
% minCliqueSize:  minimum size of cliques
% maxCliqueSize:  maximum size of cliques
% maxDegree:      the total degree of objective function
% randSeed:       random seed (need this value because this
%                 function generates polynomial randomly)
%
% <Return value>
% objPoly: objective function, which has the following fields:
% objPoly.typeCone = 1; 
% objPoly.sizeCone = 1; 
% objPoly.dimVar   --- a positive integer 
% objPoly.degree   --- a positive integer 
% objPoly.noTerms  --- a positive integer 
% objPoly.supports --- a noTerms \times nDim matrix
% objPoly.coef     --- a noTerms \times 1 vector
%
% ineqPolySys: In this case, ineqPolySys = [] because the feasible region is
%              nDim-dimensional Euclid space.
% 
% lbd: lower bound. Here, all lower bound is -1.0e+10
% ubd: upper bound. Here, all upper bound is -1.0e+10
%
% Programed by Kim sensei, and modified by H.Waki
%

%%
%% Error handling
%%
if nDim < minCliqueSize
	error('## Increase nDim more than minCliqueSize ##');
end
if mod(maxDegree,2) ~=0
	error('## maxDegree in 4th argument is even ##');
end
if nDim <= 0
	error('## nDim <= 0 ##');
elseif maxDegree < 0
	error('## maxDegree < 0 ##');
end

%%
%% Generate cliques
%%
cliqueSet = chodalStruct(nDim,minCliqueSize,maxCliqueSize,randSeed);
noOfCliques = size(cliqueSet,1);

%%
%% Initilize 'objPoly'
%%
rand('state',randSeed);
objPoly.typeCone = 1; 
objPoly.sizeCone = 1; 
objPoly.dimVar = nDim; 
objPoly.degree = maxDegree; 
objPoly.noTerms = 0; 
objPoly.supports = []; 
objPoly.coef = []; 

%%
%% If we choose a larger value than 'minEigthreshhold' as minimum
%% eigen value of V_i, the optimal sol of the given POP can be in
%% the set \{\x\in\Real^{nDim}| \max\{x_i\}<= 1\}.  
%% 
minEigthreshhold = 0;
%%
%% Generate poynomial g_i(\x) 
%%
for i=1:noOfCliques 
  subspaceIdxSet = find(cliqueSet(i,:));
  %% 1 lower degree terms
  tempSupSet = genSimplexSupport(1,nDim,maxDegree-1,subspaceIdxSet);
  [m,n] = size(tempSupSet);
  signVec = 2 * round(rand(m,1)) - 1; 
  %%
  %% Generate sparse coefficient vector
  %%
  nonZeroSW = 0; 
  while nonZeroSW == 0
    tempCoefVec = 0.001 * round(1000 * sprand(m,1,0.5)); 
    I = find(tempCoefVec); 
    nonZeroSW = ~isempty(I);
  end
  tempCoefVec = signVec .* tempCoefVec;
  objPoly.supports  = [objPoly.supports; tempSupSet(I,:)]; 
  objPoly.coef = [objPoly.coef; tempCoefVec(I)];
  %%
  %% Use 2-norm
  %%
  minEigthreshhold = minEigthreshhold + sqrt(tempCoefVec(I)'* ...
					     tempCoefVec(I)*length(I));
 
  %%
  %% Use 1-norm, in general this value is larger than one by using
  %% 2-norm.
  %%
  %minEigthreshhold = minEigthreshhold + sum(abs(tempCoefVec))*length(I);
end
%fprintf('minEigthreshhold = %g\n',minEigthreshhold); 

%%
%% Generate the highest term \v_i(\x)^T\V_i\v_i(\x)
%%                            =\Tr(\V_i \v_i(\x)\v_i(\x)^T)
%%   
for p=1:noOfCliques 
  %%
  %% Construct \v_i(\x)\v_i(\x)^T
  %%
  subspaceIdxSet = find(cliqueSet(p,:)); 	
  noMono = length(subspaceIdxSet);
  tempSupSet = sparse(zeros(noMono,nDim));
  tempSupSet(:,subspaceIdxSet) = ceil(maxDegree/2)*speye(noMono);
  %% x_i^2 (i in clique) * x_j^2 (j in clique)
  SupMat = repmat(tempSupSet,noMono,1)+kron(tempSupSet,ones(noMono,1)); 
  
  %%
  %% Construct a psd matrix V_i
  %%
  QMat = orth(2*rand(noMono,noMono)-1);		
  EMat =  minEigthreshhold*diag(rand(noMono,1)+1); 
  CMat = QMat *EMat *QMat';% coefficient matrix
  
  %%
  %% Substitute into 'objPoly'
  %%
  objPoly.supports = [objPoly.supports;SupMat];
  objPoly.coef = [objPoly.coef; reshape(CMat,noMono*noMono,1)];
end
objPoly.noTerms =  size(objPoly.supports,1);

objPoly = simplifyPolynomial(objPoly);
%%
%% Normalize 'objPoly' by dividing maximum absolute value in
%% 'objPoly'
%%
%maxVal = max(abs(objPoly.coef));
%objPoly.coef = objPoly.coef/maxVal;
%[maxVal ,maxIndex] =  max(abs(objPoly.coef));
%[minVal ,minIndex] =  min(abs(objPoly.coef));
%fprintf('maxVal = %g, minVal = %g in objPoly\n',objPoly.coef(maxIndex), objPoly.coef(minIndex));

%%
%% Set 'ineqPolySys' as empty value
%%
ineqPolySys = [];
lbd = -1.0e+10*ones(1,nDim);
ubd = 1.0e+10*ones(1,nDim);

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = genSimplexSupport(supType,nDim,r,subspaceIdxSet,order);
% supType = 0 or 1
% supspaceIdxSet \subset 1:nDim
%
% If supType == 0, then supSet consists of the elements 
% \{ v \in Z^n_+ : \sum_{i=1}^n v_i = r, v_j = 0 (j \not\in supSpaceIdxSet \}
% in the lexico graphical order. 
%
% If supType == 1 then supSet consists of the elements 
% \{ v \in Z^n_+ : \sum_{i=1}^n v_i \le r, v_j = 0 (j \not\in supSpaceIdxSet \}
% in the lexico graphical order. 
% 
if nargin < 5
  order = 'grevlex';
end
dimSubspace = length(subspaceIdxSet); 
if nDim <= 0
  error('!!! nDim <= 0 !!!');
elseif r < 0
  error('!!! r < 0 !!!');	
elseif dimSubspace == nDim
  if subspaceIdxSet(nDim) == nDim
    if supType == 0
      supSet = flatSimpSup(nDim,r);
    elseif supType == 1
      supSet = fullSimpSup(nDim,r);
    else
      error('!!! supType is neither 0 nor 1 !!!');
    end
  else
    error('!!! dimSubspace = nDim but subspaceIdxSet(nDim) not= nDim !!!');
  end
else
  supSet = restSimpSup(supType,nDim,r,subspaceIdxSet);
end

% Kojima, 02/15/2005
% [supSet,I] = Monomial_Sort(supSet, subspaceIdxSet, order);
[supSet,I] = monomialSort(supSet, subspaceIdxSet, order);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = flatSimpSup(nDim,r); 
% \{ v \in Z^n_+ : \sum_{i=1}^n v_i = r \}

if nDim == 0
  supSet = 0; 
elseif r == 0
  supSet = sparse(1,nDim);
elseif nDim == 1 
  supSet = r; 
else
  NumElem = nchoosek(nDim+r-1,r);
  supSet = sparse(NumElem,nDim);
  index = 0;
  for i=0:r
    aSupSet = flatSimpSup(nDim-1,i);
    [m,n] = size(aSupSet);
    Idx = index + [1:m];
    supSet(Idx,1) = repmat(r-i,m,1);
    supSet(Idx,2:n+1) = aSupSet;
    index = index + m;
  end
end
return; 		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = fullSimpSup(nDim,r);
% \{ v \in Z^n_+ : \sum_{i=1}^nDim v_i \leq r \}

if nDim == 1;
  supSet = [0:r]';
  supSet = sparse(supSet);
elseif r == 0
  supSet = sparse(1,nDim);
else
  NumElem = nchoosek(nDim+r,r);
  supSet = sparse(NumElem,nDim);
  index = 0;
  for i=0:r
    aSupSet = flatSimpSup(nDim,i);
    [m,n] = size(aSupSet);
    Idx = index + [1:m];
    supSet(Idx,:) = aSupSet;
    index = index + m;
  end
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function supSet = restSimpSup(supType,nDim,r,subspaceIdxSet);

dimSubspace = length(subspaceIdxSet); 
%
% error handlng
%
if (nDim < dimSubspace) | (nDim < subspaceIdxSet(dimSubspace)) 
  error('!!! nDim < dimSubspace !!!');
end
if nDim == 0;
  error('!!! nDim = 0 !!!'); 
end
if dimSubspace == 0;
  error('!!! dimSupspace = 0 !!!');
end
if nDim < subspaceIdxSet(dimSubspace)
  error('!!! nDim < subspaceIdxSet(dimSubspace) !!!'); 
end
%
% end of error handling
%

NumElem = nchoosek(dimSubspace+r, dimSubspace);
aSupSet = sparse(NumElem, dimSubspace);
if supType == 1
  %\{ v \in Z^n_+ : \sum_{i=1}^n v_i \le r, v_j = 0 (j \not\in supSpaceIdxSet) \}
  aSupSet = fullSimpSup(dimSubspace,r);
elseif supType ==0
  %\{ v \in Z^n_+ : \sum_{i=1}^n v_i = r, v_j = 0 (j \not\in supSpaceIdxSet) \}
  aSupSet = flatSimpSup(dimSubspace,r); 
else
  error('!!! You should choose 1 or 0 as supType !!!'); 
end    
m = size(aSupSet,1); 
supSet = sparse(m,nDim);
supSet(:,subspaceIdxSet) = aSupSet;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cliqueSet] = chodalStruct(nDim,minCliqueSize,maxCliqueSize,randSeed); 
%
% Generating the set of maximal cliques of a chordal graph, where 
% 1,2,\ldots,nDim is a perfect elimination ordering
%
% cliqueSet : noOfCliques \times nDim, where
% [noOfCliques,dummy] = size(cliqueSet); 
% 
% nDim = 20; 
% maxCliqueSize = 4; 
% minCliqueSize = 1; 
% randSeed = 321; 

rand('state',randSeed);
cliqueSizeVect = ceil(maxCliqueSize * rand(1,nDim)); 
cliqueSizeVect = max(cliqueSizeVect,minCliqueSize*ones(1,nDim)); 

for startIndex=1:nDim
  if startIndex==1
    [cliqueSet] = genOneClique(nDim,startIndex+1,cliqueSizeVect(1,startIndex)-1);
    cliqueSet(1,1) = 1; 
    noOfCliques = 1; 
  else 
    oneClique = sparse(zeros(1,nDim));
    oneClique(1,startIndex) = 1; 
    if startIndex < nDim
      cliqueSize = min([cliqueSizeVect(1,startIndex),nDim-startIndex+1]);
      for i=1:noOfCliques-1
	if cliqueSet(i,startIndex)==1
	  %%%
	  %%Size of one Clique may be over cliqueSize.
	  %%%
	  oneClique(1,startIndex:nDim) = max(oneClique(1,startIndex:nDim),cliqueSet(i,startIndex:nDim));
	end
      end
      I=find(oneClique);
      if length(I)>maxCliqueSize
	%fprintf('%d\n',startIndex)
	%fprintf('size of oneClique = %d\n',length(I));
	%fprintf('cliquesize        = %d\n',cliqueSize);
      end
      addSize = cliqueSize - length(find(oneClique(1,startIndex: ...
						   nDim)));
      if addSize > 0
	[addClique] = genOneClique(nDim,startIndex,addSize);
	oneClique(1,startIndex:nDim) = max(oneClique(1,startIndex:nDim),addClique(1,startIndex:nDim));
      end
      
    else
      addSize = 0;
      
    end
    
    % check maximality of oneClique ---> 
    same = 0; 
    i=1;
    sparseZeros = sparse(zeros(1,nDim-startIndex+1)); 
    while (same==0) & (i <= noOfCliques)
      diff = oneClique(1,startIndex:nDim) - cliqueSet(i,startIndex:nDim); 
      diff = max(diff,sparseZeros); 
      if length(find(diff)) == 0
	same = 1; 
      end
      i = i+1; 
    end
    % <---
    if same == 0 % oneClique is maximal compared to the cliques generated 
      if addSize >= 0
	cliqueSet = [cliqueSet; oneClique]; 
	noOfCliques = noOfCliques + 1;
      end
    end
  end
end
%sumC = sum(cliqueSet,1);
%I = find(sumC == 0);
%if length(I) > 0
%  fprintf('This case is not good because clique Set does not have all variables.\n');
%  writeCliqueSet(cliqueSet);
%end

return

function [oneClique] = genOneClique(nDim,startIndex,cliqueSize); 
% oneClique : 1 \times nDim vector
% clique = \{ i : oneClique(i) = 1 \}
if cliqueSize > nDim-startIndex+1 
  fprintf('cliqueSize = %d > nDim-startIndex+1',cliqueSize);
end
oneClique = sparse(zeros(1,nDim)); 
aVect = rand(1,nDim-startIndex+1);
[dummy,sortedIdx] = sort(aVect);
for i=1:cliqueSize
  j=sortedIdx(i) + startIndex - 1; 
  oneClique(1,j) = 1;
end
return;

function writeCliqueSet(cliqueSet);
[noOfCliques,nDim] = size(cliqueSet);
% full(cliqueSet)
fprintf('cliques\n');
for i=1:noOfCliques
  II = find(cliqueSet(i,:)); 
  cliqueSize = length(II); 
  fprintf('%2d    ',i); 
  for j=1:cliqueSize
    fprintf(' %2d',II(j));
  end
  fprintf('\n');
end
return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/example/POPformat/randomUnconst.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
