function [ Q ] = peter( Q , peterval )
% Q = peter( Q , p )
% Removes columns of Q randomly.
%
% Input:
%   Q           the array to be petered
%   p           how much to peter
%                       p<1: Q has approx. p*size(Q) elements
%                       p>=1: Q has approx p elements
%
% E.g.: peter(randn(2,30),4)
%
% Written by: tommsch, 2018

N=size(Q,2);
if peterval>=1; 
    if peterval>N; 
        return; end;
    P=peterval/N;
else
    P=peterval;
end;
B=rand(1,N);
B=B<P;
Q=Q(:,B);

end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.
