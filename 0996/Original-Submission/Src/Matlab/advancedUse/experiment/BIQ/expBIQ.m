
printFileName = 'BIQresults.csv'; 

fid = fopen(printFileName, 'a');
fprintf(fid, 'instance,optVal,method,lambda,LB,LBv,UB,time,iter(APGR),termcode,iterBP\n');
fclose(fid);

if 1
    params = {};
    params{end+1} = {'bqp100-1','bqp100-2','bqp100-3','bqp100-4','bqp100-5'};
    % params{end+1} = {'bqp100-1','bqp100-2','bqp100-3','bqp100-4','bqp100-5',...
    %    'bqp500-1','bqp500-2','bqp500-3','bqp500-4','bqp500-5'};
    params{end+1} = {'BP','sdpnalplus'};
    params{end+1} = {1e5};
    
    params = genAllComb(params);
    
    instances = params(:, 1);
    solvers = params(:, 2);
    lambdas = params(:, 3);
    
    % parfor ii = 1:length(instances)
    for ii = 1:length(instances)
        solveBIQ(printFileName,instances{ii}, solvers{ii}, lambdas{ii});
    end
end

if 1 
    params = {};
    params{end+1} = {'bqp250-1','bqp250-2','bqp250-3','bqp250-4','bqp250-5'};
    params{end+1} = {'BP'};
    params{end+1} = {1e5};
    
    params = genAllComb(params);
    
    instances = params(:, 1);
    solvers = params(:, 2);
    lambdas = params(:, 3);
    
    % parfor ii = 1:length(instances)
    for ii = 1:length(instances)
        solveBIQ(printFileName,instances{ii}, solvers{ii}, lambdas{ii});
    end
end

if 1
    params = {};
    params{end+1} = {'bqp500-1','bqp500-2','bqp500-3','bqp500-4','bqp500-5'};
    params{end+1} = {'BP'};
    params{end+1} = {1e5};
    
    params = genAllComb(params);
    
    instances = params(:, 1);
    solvers = params(:, 2);
    lambdas = params(:, 3);
    
    % parfor ii = 1:length(instances)
    for ii = 1:length(instances)
        solveBIQ(printFileName,instances{ii}, solvers{ii}, lambdas{ii});
    end
end