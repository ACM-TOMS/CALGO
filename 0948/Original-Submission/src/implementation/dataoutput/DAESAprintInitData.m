function s = DAESAprintInitData(sadata, varargin)%

printID.outfile = [];
printID.varnames = [];
printID.fcnnames = [];
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
                case {'outfile', 'varnames', 'fcnnames'}
                    printID.(name) = options.(name);
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
                printID.(str) = varargin{2*i};
            otherwise
                if isa(str, 'char')
                    error('''%s'' is a wrong argument!', str);
                else
                    error('Argument should be a string here!');
                end
        end
    end
end

if ~isempty(printID.outfile)
    try checkFileExist(printID.outfile);
    catch err, error(err.message);
    end
end

n = getSize(sadata);
if isempty(printID.varnames)
    varname = 'x';
else
    if ~isa(printID.varnames, 'cell')
        varname = printID.varnames;
    else
        m = length(printID.varnames);
        if n~=m, error('Size wrong!');
        else varname = printID.varnames;
        end
    end
end

%% Obtain Data
B = [];
if isSWP(sadata)
    IVset = getInitData(sadata)-1;
    i = find(IVset>=0);
    A = [i; IVset(i)]';
    k = 0;
    [m, n] = size(A);
    for j = 1:m
        d = 0;
        while d <= A(j,2)
            k = k+1;
            B(k,1) = i(j);
            B(k,2) = d;
            d = d+1;
        end
    end
    
    % Print what to initialize
    breakpoint = 75;
    s = '';
    if nargout==0
        s = sprintf('%s%s problem\n', s, func2str(getDAEfhandle(sadata)));
    end
    s = sprintf('%s%s\n', s, repmat('-', 1, breakpoint));
    if ~isempty(B)
        % tmp = var2string(varname, B, breakpoint);
        % s = sprintf('%s%s\n', s, repmat('-', 1, maxLen));
        s = sprintf('%sInitialization summary:\n%s\n', s, var2string(varname, B, breakpoint));
    else
        % initdataStr = sprintf('Initialization summary: (none)\n');
        % s = sprintf('%s%s\n', s, repmat('-', 1, length(initdataStr)));
        s = sprintf('%sInitialization summary: (none)\n', s);
    end
    
else
    s = sprintf('%s is structurally ill posed.\n', ...
        func2str(getDAEfhandle(sadata)));
end

if isempty(printID.outfile)
    if nargout==0, fprintf('%s', s); end
else
    fid = fopen(printID.outfile, 'w');
    fprintf(fid, '%s', s);
    fclose(fid);
end

end
