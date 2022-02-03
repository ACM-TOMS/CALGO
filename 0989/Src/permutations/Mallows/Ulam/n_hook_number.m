function [YoungTable,YoungCount,count,n_hook] = n_hook_number(partition)
% [YoungTable,YoungCount,count,n_hook]  =  n_hook_number(partition)
% n_hook_number: Given a diagram (partition), computes all the Young
% Tables, number of shapes in each Table, and the total hook number 
% These are auxiliary functions needed to compute the Ulam distance.
%
% INPUTS
% partition: The given partition specifying the shape
%
% OUTPUTS
% YoungTable: Array with the description of all the YoungTables
% YoungCount: Number of permutations in each of the YoungTables
% count:   Number of YoungTables
% n_hook:  hook number
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
% [3] J. S. Frame,  G. B. Robinson, R. M. Thrall. The hook graphs of the symmetric group. Canad. J. Math. 6, 316–325. 1954. 
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)

n = size(partition,2);
YoungTable = zeros(n); % Initialization of the Tables
YoungCount = zeros(n);

count = 0;             % The tables are filled
for j=1:n,                     
     for k=1:partition(n-j+1),          
          count = count+1;
          YoungTable(count,1:n-j+1) = 1;          
     end
end

[row,col] = find(YoungTable==1);   % Finally, YoungCounts and n_hook numbers are computed
n_hook = 1;
for i=1:size(row,1),
  YoungCount(row(i),col(i))  = sum(YoungTable(row(i),col(i):end)) + sum(YoungTable(row(i):end,col(i)))-1;
  n_hook = n_hook*YoungCount(row(i),col(i));
end  

