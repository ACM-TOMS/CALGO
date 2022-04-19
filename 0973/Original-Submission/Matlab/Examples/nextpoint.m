function [ x1 ] = nextpoint( x0 , alpha , beta)
%NEXTPOINT Determine next point in the partition of the interval [-1,1]
%   If the initial function has a singularity at alpha, then the function
%   obtained after transforming the subinterval [x1,x0] onto the interval 
%   [-1,1] will have a singularity at beta. 
%   It should hold that 1 < alpha < beta, and -1 < x0 =< 1 in order to 
%   ensure that x1 < x0. 

A = [beta,1,0;1,1,0;1,-1,1];
b = [alpha;x0;0];
X = A\b;
x1 = X(end);

end

