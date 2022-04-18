%% demo
% Demonstration of the toolbox |IPscatt|.
%
% Warning: This functions uses |setInput|, such that a directory and
% several files are created. Also variables are cleared.
%
init;
setInput;
seti.k = 50;
seti.incNb = 6;
seti.measNb = 12;
seti.nCD = 64;
seti.nOut = 2;
seti.pdaN = 50;
seti = setData(seti,0,2);
seti = setRecon(seti,0,2);
seti = recon(seti,0,2);
