function s = DAESAprintSolScheme(schemeData, varargin)

if isa(schemeData, 'SAdata')
    sadata = schemeData;
    if isSWP(sadata)
        schemeData = DAESAsolutionScheme(schemeData);
    end
elseif isa(schemeData, 'scheme')
    sadata = getSAdata(schemeData(1));
end

printSS.outfile = [];
printSS.varnames = [];
printSS.fcnnames = [];
printSS.detail = 0; 
% 0 for compact, 1 for full

n = length(varargin);
if n == 1
    options = varargin{1};
    if isstruct(options)
        names = fieldnames(options);
        for i=1:length(names)
            name = names{i};
            if ~isa(name, 'char')
                error('Argument should be a string here!');
            end
            switch name
                case {'outfile', 'varnames', 'fcnnames', 'detail'}
                    printSS.(name) = options.(name);
                otherwise
                    error('''%s'' is a wrong argument!', name);
            end
        end
    else
        error('Argument ''options'' should be a structure.\n');
    end
elseif n > 1
    if mod(n,2)==1
        error('Number of arguments is wrong!');
    end
    for i=1:n/2
        str = varargin{2*i-1};
        switch str
            case {'outfile', 'varnames', 'fcnnames'}
                printSS.(str) = varargin{2*i};
            case 'detail'
                switch varargin{2*i}
                    case 'full'
                        printSS.detail = 1;
                    case 'compact'
                    otherwise
                        if isa(varargin{2*i}, 'char')
                            error('''%s'' is a wrong argument!', varargin{2*i});
                        else
                            error('Argument should be a string here!');
                        end
                end
            otherwise
                if isa(str, 'char')
                    error('''%s'' is a wrong argument!', str);
                else
                    error('Argument should be a string here!');
                end
        end
    end
end

if ~isempty(printSS.outfile)
    try checkFileExist(printSS.outfile);
    catch err, error(err.message);
    end
end

if isSWP(sadata)
    %% Get solution data
    n = getSize(sadata);
    if isempty(printSS.varnames)
        varname = 'x';
    else
        if ~isa(printSS.varnames, 'cell')
            varname = printSS.varnames;
        else
            m = length(printSS.varnames);
            if n~=m, error('Size wrong!');
            else varname = printSS.varnames;
            end
        end
    end
    
    if isempty(printSS.fcnnames)
        fcnname = 'f';
    else
        if ~isa(printSS.fcnnames, 'cell')
            fcnname = printSS.fcnnames;
        else
            m = length(printSS.fcnnames);
            if n~=m, error('Size wrong!');
            else fcnname = printSS.fcnnames;
            end
        end
    end   
    
    breakpoint = 75;
    
    if printSS.detail==1
        s = sprintf('Solution scheme for ''%s'' problem\n', ...
        func2str(getDAEfhandle(sadata)));
    elseif printSS.detail==0
        s = sprintf('Compact solution scheme for ''%s'' problem\n', ...
        func2str(getDAEfhandle(sadata)));
    end
    % directly call printInitData function
    s = sprintf('%s%s', s, DAESAprintInitData(sadata, 'varnames', varname, ...
        'fcnnames', fcnname));
    s = sprintf('%s%s\n', s, repmat('-', 1, breakpoint));
    
    if printSS.detail==1
        printFcn = @DAESAprintScheme;
    elseif printSS.detail==0
        printFcn = @DAESAprintSchemeCompact;
    end
    
    % str = ''; maxLen = 0;
    for stage = 1:length(schemeData)
        tmp = printFcn(schemeData(stage), ...
            'varnames', varname, 'fcnnames', fcnname);
        s = sprintf('%s%s', s, tmp);        
    end   
else
    s = sprintf('%s is structurally ill posed.\n', ...
        func2str(getDAEfhandle(sadata)));
end

if isempty(printSS.outfile)
    if nargout==0, fprintf('%s', s); end
else
    fid = fopen(printSS.outfile, 'w');
    fprintf(fid, '%s', s);
    fclose(fid);
end
end