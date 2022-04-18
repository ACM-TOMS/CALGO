
printFileName = 'QAPresults.csv'; 
fid = fopen(printFileName,'a');
%fprintf(fid, 'instance,method,lambda,LB,LBv,UB,time,iter(APGR),termcode,iterBP\n');
fprintf(fid, 'instance,optVal,method,lambda,LB,LBv,UB,time,iter(APGR),termcode,iterBP\n');
fclose(fid);

if 1
    params = {};
    params{end+1} = {'chr12a', 'chr12b', 'chr12c', 'had12', 'nug12', 'rou12', 'scr12', 'tai12a', 'tai12b'};
    params{end+1} = {'BP','sdpnal','sdpnalAW'}; % {'BP'}; %%
    params{end+1} = {1e5};
    
    params = genAllComb(params);
    
    instances = params(:, 1);
    solvers = params(:, 2);
    lambdas = params(:, 3);
    
    % parfor ii = 1:length(instances)
    for ii = 1:length(instances)
        solveQAP(printFileName,instances{ii}, solvers{ii}, lambdas{ii});
    end
end

if 1
    
    %%  
    params = {};
    params{end+1} = {'chr15a', 'chr15b', 'chr15c', 'chr18a', 'chr18b',...
        'chr20a', 'chr20b', 'chr20c', 'chr22a', 'chr25a',...
        'nug20', 'nug25', 'nug30', 'tai30a', 'tai30b', 'tai35a', 'tai35b',...
        'bur26a', 'bur26b', 'bur26c', 'bur26d', 'bur26e', 'bur26f', 'bur26g', 'bur26h'};
    params{end+1} = {'BP', 'sdpnalAW'}; %% {'BP', 'sdpnal', 'sdpnalAW'};
    params{end+1} = {1e5};
    
    params = genAllComb(params);
    
    instances = params(:, 1);
    solvers = params(:, 2);
    lambdas = params(:, 3);
    
    % parfor ii = 1:length(instances)
    for ii = 1:length(instances)
        solveQAP(printFileName,instances{ii}, solvers{ii}, lambdas{ii});
    end
end

if 0
    %%
    params = {};
    params{end+1} = {'lipa40a', 'lipa40b', 'tai40a', 'tai40b', 'tho40',...
        'sko42', 'sko49', 'wil50', 'lipa50a', 'lipa50b', 'tai50a', 'tai50b'};
    params{end+1} = {'BP'}; %% {'BP', 'sdpnal', 'sdpnalAW'};
    params{end+1} = {1e5};
    
    params = genAllComb(params);
    
    instances = params(:, 1);
    solvers = params(:, 2);
    lambdas = params(:, 3);
    
    % parfor ii = 1:length(instances)
    for ii = 1:length(instances)
        solveQAP(printFileName,instances{ii}, solvers{ii}, lambdas{ii});
    end
end