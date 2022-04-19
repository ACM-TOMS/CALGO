function [ncount,FerrerShapes_Lengths] = ComputeFerrerShapes(n)
% [kdist] =  ComputeFerrerShapes(n)
% ComputeFerrerShapes: Computes all possible Ferrer shapes which are used to sample solutions 
% when the Mallows model is applied with the Ulam. Ferrer shapes of size k
% are saved in the file permus_per_shape_n_k.  The number of permutations
% per Ferrer shape for all possible shapes are saved in the file
% FerrerShapes_Lengths_n
%
% INPUTS
% n: Size of the problem (length of the permutations)
%
% OUTPUTS
% ncount: Number of FerrerShapes
% FerrerShapes_Lengths: Cell array with all FerrerShapes
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)

plist = partitions(n,[1:n],n);   % Compute all possible Ferrer Shapes (partitions of n
npart = size(plist,1);           % Number of Ferrer Shapes  

% Now we compute the hook length for each Ferrer Shape using the Hook
% Formula (i.e. n!/(prod h_\lambda(i,j)) where h_\lambda(i,j) is the number of cells in the hook)
% For efficiency reasons, we only compute (prod h_\lambda(i,j)) and leave
% the division for later,

for i=1:npart,
  [YoungTable,YoungCount,count,n_hook] = n_hook_number(plist(i,:));  
  spart(i) = n_hook;  % hook number
  scount(i) = count;  % Number of columns in the Young Table (l = n -LYS)
  AllSums(i,:) = sum(YoungTable);
  scolumn(i) = sum(AllSums(i,:)>0);
end

% Now we make the division and compute the number of permutations for each
% possible length. This information is saved in the file 
tfact = factorial(n);
spart = tfact./spart; % Here we divide the factorial by the product of h_\lambda values
square_spart = spart.^2;
ncount = zeros(1,n);
for i=1:npart,
  ncount(n-scount(i)+1) = ncount(n-scount(i)+1) + square_spart(i);   
end

% We save the Ferrer Shapes for each shape in different files
% and ordered according to the number of configurations
for i=1:n,
   pos_length = find(scount==(n-i+1));
   [val_m,pos_n] = sort(scolumn(pos_length),'ascend');
   %Ferrer_Shapes_l{i} = plist(pos_length(pos_n),:);
   Ferrer_Shapes_l = AllSums(pos_length(pos_n),1:(i));  
   eval(['save permus_per_shape_',num2str(n),'_',num2str(i),' Ferrer_Shapes_l']); 
   FerrerShapes_Lengths{i} = square_spart(pos_length(pos_n));
end
eval(['save FerrerShapes_ncounts_',num2str(n),' ncount']);
eval(['save FerrerShapes_Lengths_',num2str(n),' ncount FerrerShapes_Lengths']);

