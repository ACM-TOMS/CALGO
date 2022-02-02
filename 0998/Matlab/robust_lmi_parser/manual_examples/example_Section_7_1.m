clear;
A{1} = [0.1 0.9;0 0.1];
A{2} = [0.5 0;1 0.5];
A = rolmipvar(A,'A(\alpha)',2,1);
P = rolmipvar(2,2,'P(\alpha)',2,1);
LMIs = [[P A'*P; P*A P] >= 0];
optimize(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));
checkset(LMIs)
double(P)
