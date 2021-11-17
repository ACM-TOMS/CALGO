function s = DAESAprintSchemeCompact(schemeData, varargin)

printS.varnames = 'x';
printS.fcnnames = 'f';
n = length(varargin)/2;

if n > 0
    for i=1:n
        str = varargin{2*i-1};
        switch str
            case 'varnames'
                printS.varnames = varargin{2*i};
            case 'fcnnames'
                printS.fcnnames = varargin{2*i};
            otherwise
                if isa(str, 'char')
                    error('''%s'' is a wrong argument!', str);
                else
                    error('Argument should be a string here!');
                end
        end
    end
end

sadata = getSAdata(schemeData(1));
n = getSize(sadata);

if isempty(printS.varnames)
    varname = 'x';
else
    if ~isa(printS.varnames, 'cell')
        varname = printS.varnames;
    else
        m = length(printS.varnames);
        if n~=m, error('Size wrong!');
        else varname = printS.varnames;
        end
    end
end

if isempty(printS.fcnnames)
    fname = 'f';
else
    if ~isa(printS.fcnnames, 'cell')
        fname = printS.fcnnames;
    else
        m = length(printS.fcnnames);
        if n~=m, error('Size wrong!');
        else fname = printS.fcnnames;
        end
    end
end

breakpoint = 75;
spaceNum = 8;
initSpaceNum = 8;
adjustSpaceNum = 2;

[k, blk, blkposn, lin, iv, prov, eqn, var, solver, sadata] ...
    = getAll(schemeData);

numBlk = length(blk);

s = sprintf('k =% 3d: ', k);
delim = ' '; % used to be ',', this one seems ok
% maxLen = 0;

for i = 1:numBlk    
    F = unique((eqn{i}), 'rows');
    A = unique((var{i}), 'rows');
    if i>1, s = sprintf('%s%s', s, blanks(spaceNum)); end
    
    if isempty(F)
        tmp = var2string(varname, A, breakpoint);
        s = sprintf('%s [] : %s\n', s, tmp);
        % maxLen = max(maxLen, tmpmaxLen);
    else
        if lin(i)==1 || isempty(find(F(:,2)==0, 1)) % linear
            tmp = sprintf('[%s] : %s', var2string(fname, F), var2string(varname, A));
            tmp = breakLine(tmp, breakpoint, delim, spaceNum, initSpaceNum, adjustSpaceNum);
            s = sprintf('%s %s\n', s, tmp);
            % maxLen = max(maxLen, tmpmaxLen);
        else % can be simplified here
            % non-linear
            tmp = sprintf('[%s] : %s', var2string(fname, F), var2string(varname, A));
            tmp = breakLine(tmp, breakpoint, delim, spaceNum, initSpaceNum, adjustSpaceNum);
            s = sprintf('%s~%s\n', s, tmp);
            % maxLen = max(maxLen, tmpmaxLen);
        end
    end
end

if nargout==0
    fprintf(s)
end
end