function [objPoly, I01, Icomp, relaxOrder, params] = biqreadBP(instance,lambda,isBinary)
	if nargin==1
		lambda = 1.0e4; 
        isBinary = true;
    elseif nargin==2
        isBinary = true;
	end
	[Qmat0] = biqread([instance,'.sparse']);
    if isBinary
        [objPoly,I01,Icomp] = QmatToSlackBinQOP(Qmat0,lambda);
    else
        [objPoly,I01,Icomp] = QmatToSlackBoxQOP(Qmat0,lambda); 
    end
    relaxOrder = 1; 
    params.weightUbdObj = 0.5;
    optVal = optvalBIQ(instance); 
    params.UbdObjVal = optVal + params.weightUbdObj*abs(optVal); 
    params.UbdIX=objPoly.UbdIX; 
    params.sparseSW=0;
    params.optVal = optVal; 
end