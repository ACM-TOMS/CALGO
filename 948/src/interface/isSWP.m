function swp = isSWP(sadata)
%isSWP returns true if the DAE is structurally well-posed, and false 
%otherwise.
%swp = isSWP(sadata)
%sadata is an object returned by daeSA.
%
%See also daeSA.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
swp = DAESAisSWP(sadata);
end