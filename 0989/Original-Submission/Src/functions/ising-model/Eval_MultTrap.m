function [vals] = Eval_multTrap(vector)
% [vals] = Eval__multTrap(vector)
%  Eval__multTrap(vector): Evaluates the multi-objective trap function (two
%  objectives)
%                
% INPUT
% vector: the individual (vector) which represents a matrix of values for the spins
% all_obj: Structures with the interactions between variables pairs 
% OUTPUT
% vals: The evaluation of the uBQP function in each objective for the individual
%
% Last version 6/09/2015. Roberto Santana (roberto.santana@ehu.es)

NumbVar = size(vector,2);
vals = [0,0];
ntrapparam = 5;
for i=1:ntrapparam:NumbVar
   vals(1) = vals(1)+Trapn(vector(i:i+ntrapparam-1),ntrapparam);
end
auxvector = 1-vector;
ntrapparam = 4;
for i=1:ntrapparam:NumbVar
   vals(2) = vals(2)+Trapn(auxvector(i:i+ntrapparam-1),ntrapparam);
end






