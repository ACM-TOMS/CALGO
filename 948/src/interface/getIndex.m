function index = getIndex(sadata)
%getIndex returns the structural index of a DAE, if it is structurally well
%posed and NaN otherwise.
%index = getIndex(sadata)
%sadata is an object returned by daeSA.
%
%See also daeSA, isSWP.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
index = DAESAgetIndex(sadata);
end