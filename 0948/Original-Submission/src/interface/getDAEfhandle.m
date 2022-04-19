function h = getDAEfhandle(sadata)
%getDAEfhandle returns a handle for the function passed to daeSA
%h = getDAEfhandle(sadata)
%sadata is an object returned by daeSA.
%
%See also 
%  daeSA
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
h = DAESAgetDAEfhandle(sadata);
end
