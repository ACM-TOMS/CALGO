function n = getSize(sadata)
%getSize returns problem size, which is the value of the parameter n passed
%to daeSA.
%n = getSize(sadata)
%sadata is an object returned by daeSA.
%
%See also daeSA
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
n  = DAESAgetSize(sadata);
end