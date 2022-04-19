function DOF = getDOF(sadata)
%getDOF returns the number of degrees of freedom if the DAE is
%structurally well posed and NaN otherwise.
%DOF = getDOF(sadata) 
%sadata is an object returned by daeSA.
%
%See also daeSA, isSWP.
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
DOF = DAESAgetDOF(sadata);
end