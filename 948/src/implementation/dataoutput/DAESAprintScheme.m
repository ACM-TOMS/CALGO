function s = DAESAprintScheme(schemeData, varargin)

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
    fcnname = 'f';
else
    if ~isa(printS.fcnnames, 'cell')
        fcnname = printS.fcnnames;
    else
        m = length(printS.fcnnames);
        if n~=m, error('Size wrong!');
        else fcnname = printS.fcnnames;
        end
    end
end

space   = '    ';

breakpoint = 75;
spaceNum = 11;
initSpaceNum = 11;

[k, blk, blkposn, lin, iv, prov, eqn, var, solver] ...
    = getAll(schemeData);

numBlk = length(blk);
s = sprintf('STAGE  k = %i,', k);

if numBlk==1
    s = sprintf('%s %i block\n', s, numBlk);
else
    s = sprintf('%s %i blocks\n', s, numBlk);
end

delim = ' '; % used to be ',', this one seems ok
% maxLen = 0;

for i = 1:numBlk
    
    blkstart = blkposn{i}(1);
    blkend = blkposn{i}(end);
    s = sprintf('%s- Block %d:%d - \n', s, blkstart, blkend);
    
    % Which variables are used from previous stages. These are known
    % values.
    if ~isempty(prov{i})
        %else
        s = sprintf('%s%s%s', s, space, 'Using  ');
        A = unique((prov{i}),'rows');
        tmp = var2string(varname, A);
        tmp = breakLine(tmp, breakpoint, delim, spaceNum, initSpaceNum);
        % maxLen = max(maxLen, tmpmaxLen);
        s = sprintf('%s%s\n', s, tmp);
    end    
    
    s = sprintf('%s%s%s', s, space, 'Solve ');
    F = unique((eqn{i}), 'rows');
    A = unique((var{i}), 'rows');
    numvar = length(A(:,1));
    [numeqn, nf] = size(F);
    
    if solver(i)==-2 %isempty(F)
        s = sprintf('%s', s, 'nothing');
    else
        % if solver(i)==-1 %isempty(iv{i}) || isempty(find(F(:,2)==0, 1))
        if lin(i)==1 || isempty(find(F(:,2)==0, 1))
            s = sprintf('%slinear ', s);
        else
            s = sprintf('%snonlinear ', s);
        end
        
        numvar = length(A(:,1));
        
        if numeqn==1
            s = sprintf('%sequation', s);
        else
            s = sprintf('%s%dx%d', s, numeqn, numvar);
            s = sprintf('%s system   ', s);
        end
    end
    
    if isempty(iv{i})
        s = sprintf('%s\n', s);
    else
        %B = unique((iv{i}),'rows');
        %if A~=B
        %    error('cannot be here')
        %else % ??????
        if numeqn==0
            s = sprintf('%s %s\n', s, '(give initial value)');
        else
            if numvar==1
                s = sprintf('%s %s\n', s, '(give trial value)');
            else
                s = sprintf('%s %s\n', s, '(give trial values)');
            end
        end
        %end
    end
    
    if numeqn > 0 % && numeqn+numvar < 20
        tmp = sprintf('    0 = %s  for  %s', var2string(fcnname, F), var2string(varname,A));
        tmp = breakLine(tmp, breakpoint, delim, 8, 0);
        s = sprintf('%s%s\n', s, tmp);
        % maxLen = max(maxLen, tmpmaxLen);
    else
        tmp = sprintf('    for    %s', var2string(varname,A));
        tmp = breakLine(tmp, breakpoint, delim, 8, 0);
        s = sprintf('%s%s\n', s, tmp);
        % maxLen = max(maxLen, tmpmaxLen);
    end
end

if nargout==0
    fprintf(s)
end
end
