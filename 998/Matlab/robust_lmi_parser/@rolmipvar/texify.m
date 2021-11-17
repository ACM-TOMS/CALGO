function varargout = texify(varargin)
% Returns the LaTeX math code of the input polynomial
%
% out = texify(poly,option,var1,...,varN,simplexnames) returns the LaTeX
% math code of the polynomial given in poly. The inputs var1, ..., varN
% correspond to the rolmipvar variables whose labels are depicted in expr.
% If some term is not informed, then it is supposed to be a parameter
% independent variable. The optional input simplexnames is a cell array of
% strings containing the simplex names. If, instead of a cell array, just a
% string is informed, then all the unit simplexes are grouped into one
% multi-simplex variable. The options may be:
%
% 'explicit': It is shown the dependence of the variables on the simplexes
% in the general form, not showing each monomial.
%
% 'implicit': The matrices are presented without their dependence on the
% simplexes.
% 
% 'polynomial': The general polynomial, showing each monomial, is
% presented.

if (nargin < 3)
    error('Input error. Type ''help texify'' for more details');
end

poly = varargin{1};
option = varargin{2};
for cont=3:length(varargin)-1
    var{cont-2} = varargin{cont};
end
multisimplex = 0;
if (iscell(varargin{nargin}))
    simplexnames = varargin{nargin};
    numsimplexes = length(simplexnames);
else
    if (isstr(varargin{nargin})) %Multi-simplex input
        multisimplex = 1;
        simplexnames = varargin{nargin};
    else
        simplexnames = [];
        var{nargin-2} = varargin{nargin};
        numsimplexes = 0;
        for cont = 1:length(var)
            numsimplexes = max([numsimplexes length(var{cont}.data(1).exponent)]);
        end
        for cont=1:numsimplexes
            simplexnames{cont} = char('a'+cont-1);
        end
    end
end
expr = poly.label;

if (strcmp(option,'explicit'))
    mode = 1;
elseif (strcmp(option,'polynomial'))
    mode = 2;
elseif (strcmp(option,'implicit'))
    mode = 3;
else
    error('Option not recognized. Type ''help texify'' for more details');
end

%Taking off the space chars
ind = find(expr==' ');
expr(ind) = [];

out = [];
if ((mode == 1) || (mode == 3))
    cont = 1;
    if (expr(1) == '[') %Matrix of matrices
        out = ['\begin{bmatrix} '];
        cont = 3;
    end
    
        
    while (cont <= length(expr))
        varnow = [];
        while ((cont <= length(expr)) && (expr(cont)~='+') && (expr(cont)~='*') && (expr(cont)~='-') && (expr(cont)~=char(39)) && (expr(cont)~='(') && (expr(cont)~=')') && (expr(cont)~=',') && (expr(cont)~=']'))
            varnow = [varnow expr(cont)];
            cont = cont + 1;
        end
        contvar = 1;
        while ((contvar <= length(var)) && (~strcmp(var{contvar}.label,varnow)))
            contvar = contvar + 1;
        end
        if (contvar <= length(var)) %If found the variable on varnow
            out = [out var{contvar}.label];
            if (mode == 1)
                if (sum(var{contvar}.vertices) > 0)
                    out = [out '('];
                    first = 1;
                    if (multisimplex)
                        out = [out simplexnames];
                    else
                        for contsimplex = 1:length(var{contvar}.vertices)
                            if (var{contvar}.vertices(contsimplex) > 0)
                                if (~first)
                                    out = [out ','];
                                end
                                first = 0;
                                out = [out simplexnames{contsimplex}];
                            end
                        end
                    end
                    out = [out ')'];
                end
            end
        else
            out = [out varnow];
        end
        
        if (cont <= length(expr)) %Operation
            if (expr(cont) == ',') %End of cell
                out = [out ' & '];
            elseif (expr(cont) == ']') %End of row
                cont = cont + 1;
                if (expr(cont) == ';') %Continue row below
                    out = [out ' \\ '];
                    cont = cont + 1;
                elseif (expr(cont) == ']') %End of matrix
                    out = [out ' \end{bmatrix}'];
                end
            elseif (expr(cont) ~= '*')
                out = [out expr(cont)];
            end
        end
        cont = cont + 1;
    end
end

if (mode == 2)
    for cont=1:length(poly.opcode)
        strmain = [];
        if (cont > 1)
            strmain = [' + '];
        end
        for contsimplex=1:length(poly.data(cont).exponent)
            expon = poly.data(cont).exponent{contsimplex};
            for contexp=1:length(expon)
                if (expon(contexp) == 1)
                    if (multisimplex)
                        strmain = [strmain simplexnames '(' num2str(contsimplex) ')' '_{' num2str(contexp) '}'];
                    else
                        strmain = [strmain simplexnames{contsimplex} '_{' num2str(contexp) '}'];
                    end
                elseif (expon(contexp) > 1)
                    if (multisimplex)
                        strmain = [strmain simplexnames '(' num2str(contsimplex) ')' '_{' num2str(contexp) '}^{' num2str(expon(contexp)) '}'];
                    else
                        strmain = [strmain simplexnames{contsimplex} '_{' num2str(contexp) '}^{' num2str(expon(contexp)) '}'];
                    end
                end
            end
        end
        
        
        exprop = poly.opcode{cont};
        contop = 1;
        matrices = 0;
        if (exprop(1) == '[') %Matrix of matrices
            contop = 3;
            matrices = 1;
            strmain = [strmain '\left(\begin{bmatrix} '];
        else
            strmain = [strmain '('];
        end
        while (contop <= length(exprop))
            varnow = [];
            while (exprop(contop) ~= '#')
                varnow = [varnow exprop(contop)];
                contop = contop + 1;
            end
            
            varnow(find(varnow=='*')) = [];
            strmain = [strmain varnow];
            contop = contop + 2;
            if ((exprop(contop-1) == 'V') || (exprop(contop-1) == 'C'))
                %The only codes that can depend on parameters
                strmain = [strmain '_{'];
                while ((contop <= length(exprop)) && (~isnan(str2double(exprop(contop)))) )
                    %While is a number
                    strmain = [strmain exprop(contop)];
                    contop = contop + 1;
                end
                strmain = [strmain '}'];
            end

            while ((contop <= length(exprop)) && ((exprop(contop)=='+') || (exprop(contop)=='*') || (exprop(contop)=='-') || (exprop(contop)==char(39)) || (exprop(contop)=='(') || (exprop(contop)==')')))
                if (exprop(contop) ~= '*') %Don't show the * on LaTex!
                    strmain = [strmain exprop(contop)];
                end
                contop = contop + 1; %Find the operation
            end            
            
            if (contop <= length(exprop))
                if (exprop(contop) == ',') %End of cell
                    strmain = [strmain ' & '];
                    contop = contop + 1;
                elseif (exprop(contop) == ']') %End of row
                    contop = contop + 1;
                    if (exprop(contop) == ';') %Continue row below
                        strmain = [strmain ' \\ '];
                        contop = contop + 1;
                    elseif (exprop(contop) == ']') %End of matrix
                        strmain = [strmain ' \end{bmatrix}'];
                    end
                    contop = contop + 1;
                end
            end
            
        end
        if (matrices)
            strmain = [strmain '\right)' char(10)];
        else
            strmain = [strmain ')' char(10)];
        end
        out = [out strmain];
    end
end
             
varargout{1} = out;
    


return
