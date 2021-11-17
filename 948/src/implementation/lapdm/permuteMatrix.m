function [A, p, q, r, s, cc, rr] = permuteMatrix(A)

Asp = sigmaToSparse(A);
[p, q, r, s, cc, rr] = dmperm(Asp); 
A = A(p, q); 

end