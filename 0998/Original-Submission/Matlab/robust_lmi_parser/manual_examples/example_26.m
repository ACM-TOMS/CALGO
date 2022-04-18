clear;

A{1} = {[1 0],eye(2)};
A{2} = {[0 1],2*eye(2)};

polyA = rolmipvar(A,'A',2,1);
polyP = rolmipvar(2,2,'P','sym',2,1);
cond = [polyP polyA'*polyP; polyP*polyA polyP];
LMIs = [cond > 0]
d = 3;
condpolya = polya(cond,d);
LMIspolya = [condpolya > 0]
