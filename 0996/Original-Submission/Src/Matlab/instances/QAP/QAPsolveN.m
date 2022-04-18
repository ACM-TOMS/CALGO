function [sol,info]=QAPsolveN(paramBP,instance)

filename = [instance '.dat'];
optVal=optvalQAP(instance); 
[objPoly,Icomp,I0,penaltyPoly,A,mb] = qapreadN(filename);

nn = size(objPoly.supports, 2);

penaltyPoly.coef = ...
    paramBP.lambda*(max(1e-8,norm(objPoly.coef))/max(1e-8, norm(penaltyPoly.coef)))*penaltyPoly.coef;

objPoly = addPoly(objPoly,penaltyPoly); 

% polySys = [];
switch paramBP.method
    case 1
        % complementarity is in polyCone
        % box
        % Icomp2 = Icomp; % complementarity is in polyCone
        I0 = ~I0; % box
    case 2
        % no complementarity
        % binary
        Icomp = []; % no complementarity
        % I012 = I01; % binary 
    case 3
        % complementarity is in polyCone
        % binary
        % Icomp2 = Icomp; % complementarity is in polyCone
        % I012 = I01; % binary
end

paramBP.UbdObjVal = optVal + paramBP.weightUbdObj*abs(optVal); 
paramBP.UbdIX=sqrt(nn)+1;
paramBP.sparseSW=0;

[sol,info] = BPall(objPoly,I0,Icomp,paramBP);

end

