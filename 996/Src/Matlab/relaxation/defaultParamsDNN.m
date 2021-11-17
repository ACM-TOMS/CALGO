function params = defaultParamsDNN()
%DEFAULTPARAMSDNN   Default parameters for BBPOPtoCOP

    params.printyesDNN = 2; 
    
%   This function returns default parameters.

    params.sparseSW = true;
    params.DNNSW = true;
    params.relaxOrder = 1;

%     params.COPrep = 'basis';
%     params.lambda = 1e3;
%     params.maxtime = 1e5;
%     params.solver = 'BP';
    
    params.LagDNN = false; % strengthen relaxation
    
    params.findLongChain = true;
    params.findArborescence = false;
end

