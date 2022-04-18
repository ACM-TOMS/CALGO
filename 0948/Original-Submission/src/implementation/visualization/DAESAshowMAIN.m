function DAESAshowMAIN(sadata, varargin)

if ~isa(sadata, 'SAdata')
    error('The first argument should be an SAdata object.\n');
end

showFcn = @DAESAshowStruct;
sm = [{0}, {[]}];

n = length(varargin);
if n == 1 && isstruct(varargin{1})
    options = varargin{1};
    names = fieldnames(options);
    smSet = 0; % in case we have both 'submat' and 'blocksubmat'
    for i=1:length(names)
        name = names{i};
        switch name
            case 'disptype'
                switch options.disptype
                    case 'original'
                    case 'blocks'
                        showFcn = @DAESAshowStructBlks;
                    case 'fineblocks'
                        showFcn = @DAESAshowStructFineBlks;
                    otherwise
                        error('''disptype'' should be ''original'', ''blocks'' or ''fineblocks''.')
                end
            case 'submat'
                if smSet==1, error('Cannot have both submat and blocksubmat!'); end
                sm = [{1}, {options.submat}];
                smSet = 1;
            case 'blocksubmat'
                if smSet==1, error('Cannot have both submat and blocksubmat!'); end
                if isempty(options.disptype) || strcmp(options.disptype, 'original')
                    error('For original display, it cannot display blocks!');
                end
                sm = [{2}, {options.blocksubmat}];
                smSet = 1;
            otherwise
                error('''%s'' is a wrong argument!', name);
        end
    end
    
else
    if mod(n,2)==1
        error('Number of arguments is wrong or the argument is not a structure!');
    end
    smSet = 0; % in case we have both 'submat' and 'blocksubmat'
    for i=1:n/2
        optArgName = varargin{2*i-1};
        optArgVal = varargin{2*i};
        switch optArgName
            case 'disptype'
                switch optArgVal
                    case 'original'
                    case 'blocks'
                        showFcn = @DAESAshowStructBlks;
                    case 'fineblocks'
                        showFcn = @DAESAshowStructFineBlks;
                    otherwise
                        error('''disptype'' should be ''original'', ''blocks'' or ''fineblocks''.')
                end
            case 'submat'
                if smSet==1, error('Redundant arguments!'); end
                sm = [{1}, {optArgVal}];
                smSet = 1;
            case 'blocksubmat'
                if smSet==1, error('Redundant arguments!'); end
                sm = [{2}, {optArgVal}];
                smSet = 1;
            otherwise
                if isa(optArgName, 'char')
                    error('''%s'' is a wrong argument!', optArgName);
                else
                    error('Argument should be a string here!');
                end
        end
    end
    if strcmp(func2str(showFcn), 'DAESAshowStruct')
        if sm{1}==1
            sm = {sm{2}}; showFcn = @DAESAshowSM;
        elseif sm{1}==2
            error('For original display, it cannot display blocks!');
        end
    end
end

if smSet==0
    showFcn(sadata);
elseif smSet==1
    showFcn(sadata, sm{:});
end
end