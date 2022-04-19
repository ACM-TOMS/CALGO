% GETDEGREE   GETDEGREE(OBJ) provides the apparent degree of the
% polynomial (the true degree can be lower)
function n = getDegree(obj)
    n = length(obj.BernsCoeff)-1;
end


