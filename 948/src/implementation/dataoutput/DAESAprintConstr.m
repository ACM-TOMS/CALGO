function s = DAESAprintConstr(sadata, varargin)

printC.outfile = [];
printC.fcnnames = [];
printC.varnames = [];

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
                    printC.(name) = options.(name);
                otherwise
                    error('''%s'' is a wrong argument!', name);
            end
        end
    else
        error('Argument ''options'' should be a structure.');
    end
elseif n > 1
    if mod(n,2)==1
        error('Number of arguments is wrong!');
    end
    for i=1:n/2
        str = varargin{2*i-1};
        switch str
            case {'outfile', 'varnames', 'fcnnames'}
                printC.(str) = varargin{2*i};
            otherwise
                if isa(str, 'char')
                    error('''%s'' is a wrong argument!', str);
                else
                    error('Argument should be a string here!');
                end
        end
    end
end

if ~isempty(printC.outfile)
    try checkFileExist(printC.outfile);
    catch err, error(err.message);
    end
end

n = getSize(sadata);

if isempty(printC.fcnnames)
    fcnname = 'f';
else
    if ~isa(printC.fcnnames, 'cell')
        fcnname = printC.fcnnames;
    else
        m = length(printC.fcnnames);
        if n~=m, error('Size wrong!');
        else fcnname = printC.fcnnames;
        end
    end
end

%% Obtain Data
B = [];
if isSWP(sadata)
    constr = getConstr(sadata)-1;
    i = find(constr>=0);
    A = [i; constr(i)]';
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
    
    breakpoint = 75;
    % Print the constraints
    s = sprintf('%s problem\n', func2str(getDAEfhandle(sadata)));
    s = sprintf('%s%s\n', s, repmat('-', 1, breakpoint));
    if ~isempty(B)
        % tmp = var2string(constrname, B, breakpoint);
        % s = sprintf('%s%s\n', s, repmat('-', 1, breakpoint));
        s = sprintf('%sConstraints:\n%s\n', s, var2string(fcnname, B, breakpoint));
    else
        % constrStr = sprintf('Constraints: (none)\n');
        s = sprintf('%sConstraints: (none)\n', s);
        % s = sprintf('%s%s', s, constrStr);
    end
else
    s = sprintf('%s is structurally ill posed.\n', ...
        func2str(getDAEfhandle(sadata)));
end

if isempty(printC.outfile)
    if nargout==0, fprintf('%s', s); end
else
    fid = fopen(printC.outfile, 'w');
    fprintf(fid, '%s', s);
    fclose(fid);
end
end
