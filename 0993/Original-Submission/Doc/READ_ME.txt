ckronx is a self contained Matlab function to perform the operation
kron(A1,A2, ... ,Ad)*B
where Ai is an mi x ni matrix and B is an n1n2...nd x q matrix. The result is an m1m2...md x q matrix. 

type 
  help ckronx
at the Matlab command line for detailed documentation on usage.

Example:
   m=[50 50]; n=[30 60]; 
   d=length(m);
   A=cell(1,d);
   for i=1:d, A{i}=randn(m(i),n(i)); end
   B=randn(prod(m),5);
   options=struct('transpose',ones(1,d));
   C=ckronx(A,B,options);