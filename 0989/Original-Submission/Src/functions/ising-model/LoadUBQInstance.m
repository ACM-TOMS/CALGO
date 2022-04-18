function [all_obj] = LoadUBQInstance(fname,n,nobj)
% [all_obj] = LoadUBQInstance(fname,nobj)
% LoadIsing: Reads  the file of a multi-objective unconstrained UBQ problem
% INPUT
% n: number of variables
% nobj: number of objectives
% OUTPUT
% 
% all_obj: all_obj{i} contains all the interactions for object i.
%
% Last version 6/09/2015. Roberto Santana (roberto.santana@ehu.es)



fid = fopen(fname,'r');
seed = fscanf(fid,'%d',1);
nedges = zeros(1,nobj);
for j=1:nobj,
  all_obj{j} = [];      
  A = fscanf(fid,'%d',2);
  nedges(j) = A(2);
  for i=1:nedges(j),
    all_obj{j}(i,:) = fscanf(fid,'%d',3);    
  end     
end
fclose(fid);

