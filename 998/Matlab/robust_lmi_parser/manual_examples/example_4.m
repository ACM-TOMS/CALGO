clear;

ma{1} = {[1 0],-1};
ma{2} = {[0 1],1};
a = rolmipvar(ma,'a',2,1);
mb{1} = {[0 0],[1 0],-3};
mb{2} = {[0 0],[0 1],2};
b = rolmipvar(mb,'b',[0 2],[0 1]);


% %Alternative way using ''fork''
% ma{1} = {[1 0],-1};
% ma{2} = {[0 1],1};
% a = rolmipvar(ma,'a',2,1);
% mb{1} = {[1 0],-3};
% mb{2} = {[0 1],2};
% b = rolmipvar(mb,'b',2,1);
% b = fork(b,'b',1,2)
