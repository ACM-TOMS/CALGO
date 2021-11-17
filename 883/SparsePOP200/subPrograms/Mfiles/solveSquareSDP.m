function [x,y,infeasibleSW] = solveSquareSDP(A,b,c,K); 
if ~isempty(K.s)
    fprintf('## The coeeficient matrix A is not legitimate! ##\n'); 
    exit;
end
infeasibleSW = -1; 
x = A\b; y = A'\c; 
pointer = K.f; 
if K.l > 0 
    i = pointer;
    while (i < pointer+K.l) & (infeasibleSW == -1) 
    	i=i+1;
        if x(i) <= -1.0e-6
            infeasibleSW = 1;
        end
    end
    pointer = pointer + K.l;
end

return
