function poly = rolmipvar(varargin)
% Procedure to compose the polynomial structs to use the parser.
%
% [poly] = rolmipvar(M,label,vertices,degree) returns the polynomial structure
% for a variable set M, with a given label and given the number of 
% vertices and the degree (degree = 0 if constant). The variables
% M1, M2, ... , MN can be concatenated in M (M = [M1 M2 ... MN]) or
% disposed in cells, either inputting the coefficients without informing 
% the exponents of each monomial (M{1} = M1, ... , M{N} = MN) (in this
% case, the user is supposed to know the order of the exponents, which is
% the same order of the output of the function
% generate_homogeneous_exponents), either informing the exponents of each
% monomial (M{1} = {[1 0 0 ... 0 0], M1}, ... M{N} = {[0 0 0 ... 0 1],MN}).
%
% poly = rolmipvar(M,label,'scalar') returns a scalar. Note that 
% M has to have dimension 1x1.
%
% poly = rolmipvar(rows,cols,label,vertices,degree) internally 
% defines a set of rows X cols matrix variables M (symmetric if square,
% full otherwise) with a given label and given the number of vertices and
% the degree.
%
% poly = rolmipvar(rows,cols,label,parametr,vertices,degree)
% defines a variable: symmetric if parametr = 'symmetric'; full if parametr
% = 'full'; symmetric toeplitz if parametr = 'toeplitz'; symmetric hankel 
% if parametr = 'hankel'; and skew-symmetric if parametr = 'skew'.
%
% poly = rolmipvar(label,'scalar') internally defines a scalar
% variable.
%
% Special cases (important to use with the lmifiles command): 
%    rolmipvar(eye(n),label) returns the identity matrix with a label
%    rolmipvar(zeros(n,m),label) returns the zero matrix with a label
%
% The command rolmipvar may also be used to define matrices polynomially
% dependent on parameters with known bounds.
%
% [poly] = rolmipvar(M,label,bounds) returns the polynomial structure,
% already converted to simplex representation (each parameter generates a
% different simplex), for a variable set M, with a given label and the
% known bounds. In this case, the input M must be a cell array variable, each
% cell containing a monomial. For instance, if the polinomial is equal to
% A(p) = A0 + A1*p1 + A2*p2^2, matrix M is defined as
% M{1} = {[0 0],A0}; M{2} = {[1 0],A1}; M{2} = {[0 2],A2};
%
% [poly] = rolmipvar(rows,cols,label,parametr,polmask,bounds) defines
% matrix variables of dimensions rows X cols, with structure given by
% parametr (is "symmetric" if not informed), whose dependence on the 
% parameters (with known bounds) are defined by polmask. The input polmask 
% is a cell array variable. For instance, if the variable to be defined is 
% given by V(p) = V0 + V1*p1 + V2*p2^2, then
% polmask = {[0 0],[1 0],[0 2]}.

scalar = false;
polynomial = false;
if (nargin == 1)
  M = varargin{1};
  if(isstruct(M))
    poly.label = M.label;
    poly.vertices = M.vertices;
    poly.opcode = M.opcode;
    poly.data = M.data;
    if (isfield(M,'bounds'))
        poly.bounds = M.bounds;
    else
        poly.bounds = [];
    end
    poly = class(poly, 'rolmipvar');
    return;
  end
elseif (nargin == 2)
    M = varargin{1};
    label = varargin{2};
    vertices = 0;
    degree = 0;
    if (strcmp(label,'scalar'))
        label = M;
        scalar = true;
    end
elseif (nargin == 3) 
    M = varargin{1};
    if (isstr(varargin{3}))
        if (length(M) > 1)
            error('The first input argument should be a scalar');
            poly = [];
            return
        end
        label = varargin{2};
        vertices = 0;
        degree = 0;
        if (strcmp(varargin{3},'scalar'))
            scalar = true;
        else
            error('Invalid third argument');
            poly = [];
            return
        end
    else %Polynomial parameter representation
        polynomial = true;
        label = varargin{2};
        bounds = varargin{3};
    end
elseif (nargin == 4)
    M = varargin{1};
    label = varargin{2};
    vertices = varargin{3};
    degree = varargin{4};
elseif (nargin == 5)
    if (~iscell(varargin{4}))
        rows = varargin{1};
        cols = varargin{2};
        label = varargin{3};
        vertices = varargin{4};
        degree = varargin{5};
        parametr = 'symmetric';
        parametrextra = 'real';
    else
        rows = varargin{1};
        cols = varargin{2};
        label = varargin{3};
        polmask = varargin{4};
        bounds = varargin{5};
        parametr = 'symmetric';
        parametrextra = 'real';
    end
elseif (nargin == 6)
    if (~iscell(varargin{5}))
        rows = varargin{1};
        cols = varargin{2};
        label = varargin{3};
        parametr = varargin{4};
        parametrextra = 'real';
        vertices = varargin{5};
        degree = varargin{6};
    else %Polynomial parameter representation
        polynomial = true;
        rows = varargin{1};
        cols = varargin{2};
        label = varargin{3};
        parametr = varargin{4};
        parametrextra = 'real';
        polmask = varargin{5};
        bounds = varargin{6};
    end
elseif (nargin == 7)
    if (~iscell(varargin{5}))
        rows = varargin{1};
        cols = varargin{2};
        label = varargin{3};
        parametr = varargin{4};
        parametrextra = varargin{5};
        vertices = varargin{6};
        degree = varargin{7};
    else %Polynomial parameter representation
        polynomial = true;
        rows = varargin{1};
        cols = varargin{2};
        label = varargin{3};
        parametr = varargin{4};
        parametrextra = varargin{5};
        polmask = varargin{6};
        bounds = varargin{7};
    end
end

if (polynomial)
    if (nargin ~= 3)
        %Define the variable
        for cont=1:length(polmask)
            M{cont}{1} = polmask{cont};
            M{cont}{2} = sdpvar(rows,cols,parametr,parametrextra);
        end
    end
    numparam = length(M{1}{1});
    
    constant = [];
    polyend = [];
    for cont=1:length(M)
        if (sum(M{cont}{1}) == 0)
            constant = M{cont};
        else
            ind = find(M{cont}{1});
            poly = [];    
            for contmon=1:length(ind)
                %Determine the coefficients
                Maux = [];
                deg = M{cont}{1}(ind(contmon));
                coefs = [bounds(ind(contmon),1) bounds(ind(contmon),2)];
                for contcoef = 2:deg
                    coefs = conv(coefs,[bounds(ind(contmon),1) bounds(ind(contmon),2)]);
                end
                
                expon = deg;
                for contcoef = 1:length(coefs)
                    Maux{contcoef}{ind(contmon)} = [expon deg-expon];
                    Maux{contcoef}{numparam + 1} = coefs(contcoef);
                    expon = expon - 1;
                    
                    %Zero the other simplexes
                    for contaux = 1:numparam
                        if (contaux ~= ind(contmon))
                            Maux{contcoef}{contaux} = 0;
                        end
                    end
                end
                
                vert = zeros(1,numparam);
                vert(ind(contmon)) = 2;
                degr = zeros(1,numparam);
                degr(ind(contmon)) = deg;
                polyaux = rolmipvar(Maux,'1',vert,degr);
                if (isempty(poly))
                    poly = polyaux;
                else
                    poly = poly*polyaux;
                end
            end
            poly = poly*M{cont}{2};
            if (isempty(polyend))
                polyend = poly;
            else
                polyend = polyend + poly;
            end
        end
    end
    poly = polyend;
    
    Maux = [];
    if (~isempty(constant))
%         for contaux = 1:numparam
%             Maux{1}{contaux} = 0;
%         end
%         Maux{1}{numparam+1} = constant{2};
        Maux = constant{2};
        polyaux = rolmipvar(Maux,'0',0,0);
        poly = poly + polyaux;
    end
    
    poly.bounds = bounds;
    poly.label = label;
else
    numtotal = 1; %numcoefs total
    if (nargin >= 5)
        %Define the variable
        conttotal = 1;
        for (contsimplex = 1:length(vertices))
            if (vertices(contsimplex) > 0)
                numcoefs = (factorial(vertices(contsimplex) + ...
                    degree(contsimplex) - 1))/(factorial(degree(contsimplex)) ...
                    *factorial(vertices(contsimplex)-1));
            else
                numcoefs = 1;
            end
            if (degree(contsimplex) == 0)
                numcoefs = 1;
            end
            numtotal = round(numtotal*numcoefs);
        end
        for cont=1:numtotal
            M{conttotal} = sdpvar(rows,cols,parametr,parametrextra);
            conttotal = conttotal + 1;
        end
    end
    
    if ((nargin == 2) && (scalar))
        M = sdpvar(1,1);
    end
    
    
    poly.label = label;
    poly.vertices = vertices;
    
    if (~iscell(M))
        if (length(vertices) > 1)
            error('In the multi-simplex case, input each monomial as {deg1,deg2,...,degN,A}.');
            poly = [];
            return
        else
            if (degree > 0) %Not a constant
                numcoefs = (factorial(vertices + degree - 1))/(factorial(degree)*factorial(vertices-1));
                [order,tam] = size(M);
                tam = tam/numcoefs;
                [exponents] = generate_homogeneous_exponents(vertices,degree);
                for cont=1:size(exponents,1)
                    poly.data(cont).exponent = {exponents(cont,:)};
                    poly.data(cont).value = M(:,(cont-1)*tam+1:cont*tam);
                end
            else
                % Alexandre Felipe 2015-Mar-27
                %    Even if we create a degree zero polynomial it must
                %    have the exponents of all variables of the simplex.
                %    {zeros(1, vertices)} instead of {0}
                poly.data(1).exponent = {zeros(1, max([vertices 1]))};
                poly.data(1).value = M;
            end
        end
        
    else %If M is a cell array
        if (~iscell(M{1})) %If the degree of each monomial is not specified
            if ((length(vertices) > 1) && (nargin < 5))
                error('In the multi-simplex case, input each monomial as {deg1,deg2,...,degN,A}.');
                poly = [];
                return
            else
                if (nargin < 5) %The matrices are given in the input data
                    if (degree > 0)
                        numcoefs = (factorial(vertices + degree - 1))/(factorial(degree)*factorial(vertices-1));
                        if (length(M) ~= numcoefs)
                            error('Invalid number of monomials entered');
                            poly = [];
                            return
                        else
                            [exponents] = generate_homogeneous_exponents(vertices,degree);
                            for cont=1:size(exponents,1)
                                poly.data(cont).exponent = {exponents(cont,:)};
                                poly.data(cont).value = M{cont};
                            end
                        end
                    else
                        poly.data(1).exponent{1} = zeros(1,max([vertices 1]));%0; BATEMA
                        poly.data(1).value = M{1};
                    end
                else %A variable is being generated
                    numcoefs = 0;
                    for contsimplex=1:length(degree)
                        if (degree(contsimplex) > 0) %Not a constant
                            [exponents{contsimplex}] = generate_homogeneous_exponents(vertices(contsimplex),degree(contsimplex));
                            numcoefs = numcoefs + (factorial(vertices(contsimplex) + degree(contsimplex) - 1))/(factorial(degree(contsimplex))*factorial(vertices(contsimplex)-1));
                        else %degree(contsimplex) == 0
                            exponents{contsimplex} = zeros(1,max([vertices(contsimplex) 1])); %exponents{contsimplex} = 0; BATEMA
                            numcoefs = numcoefs + 1;
                        end
                        exponent{contsimplex} = exponents{contsimplex}(1,:);
                    end
                    
                    if (numtotal > 1)
                        numcoefs = numtotal;
                    end
                    
                    contexp = ones(1,length(vertices)+1);
                    if ((nargin < 5) && (length(M) ~= numcoefs))
                        error('Invalid number of monomials entered');
                        poly = [];
                        return
                    else
                        for cont=1:length(M)
                            poly.data(cont).exponent = exponent;
                            poly.data(cont).value = M{cont};
                            
                            contexp(1) = contexp(1) + 1;
                            if (size(exponents{1},1) >= contexp(1))
                                exponent{1} = exponents{1}(contexp(1),:);
                            end
                            aux = 1;
                            while ((aux <= length(vertices)) && (contexp(aux) > size(exponents{aux},1)))
                                contexp(aux) = 1;
                                contexp(aux+1) = contexp(aux+1) + 1;
                                
                                exponent{aux} = exponents{aux}(1,:);
                                if ((aux < length(exponents)) && (size(exponents{aux+1},1) >= contexp(aux+1)))
                                    exponent{aux+1} = exponents{aux+1}(contexp(aux+1),:);
                                end
                                
                                aux = aux + 1;
                            end
                        end
                    end
                end
            end
        else  %If the degree of each monomial is specified
            %Create the hash table
            jump = 1;
            numelem = 1;
            for contsimplex=1:length(vertices)
                [exponents{contsimplex}] = generate_homogeneous_exponents(vertices(contsimplex),degree(contsimplex));
                if (contsimplex > 1)
                    jump(contsimplex) = numelem;
                end
                numelem = numelem*size(exponents{contsimplex},1);
            end
            
            for cont=1:length(M)
                if (length(vertices) ~= length(M{cont})-1)
                    error('The definition of the matrices does not correspond to the informed number of simplexes');
                else
                    for contsimplex=1:(length(M{cont})-1)
                        if (sum(M{cont}{contsimplex}) == degree(contsimplex))
                            exponent{contsimplex} = M{cont}{contsimplex};
                        else
                            error('At least one exponent of the polynomial does not have the chosen degree');
                            poly = [];
                            return
                        end
                    end
                end
                
                indresul = gethash(exponent,exponents,jump);
                poly.data(indresul).exponent = exponent;
                poly.data(indresul).value = M{cont}{length(M{cont})};
            end
            
            %Insert the zero-monomials
            
            contexp = ones(1,length(vertices)+1);
            for conttotal=1:numelem
                exponent = [];
                for contsimplex=1:(length(M{cont})-1)
                    exponent{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);
                end
                
                indresul = gethash(exponent,exponents,jump);
                if ((indresul > length(poly.data)) || (isempty(poly.data(indresul).value)))
                    poly.data(indresul).exponent = exponent;
                    poly.data(indresul).value = zeros(size(M{1}{length(M{1})}));
                end
                contexp(1) = contexp(1) + 1;
                aux = 1;
                while ((aux <= length(vertices)) && (contexp(aux) > size(exponents{aux},1)))
                    contexp(aux) = 1;
                    contexp(aux+1) = contexp(aux+1) + 1;
                    aux = aux + 1;
                end
            end
        end
    end
    
    
    %Determine the field opcode: if the polynomials are variables (V),
    %constants (C), predefined identity or zeros matrices (P), scalar
    %variable (S) or scalar constant (K)
    %The type matrix composition (M) may appear through the construct_lmi_terms
    if (isa(poly.data(1).value,'sdpvar'))
    %if (sum(sum(isnan(double(poly.data(1).value)))) > 0)
        %if (isnan(double(poly.data(1).value(1,1))))
        if (scalar)
            poly.opcode{1} = strcat(poly.label,'#S1');
        else
            for cont=1:length(poly.data)
                poly.opcode{cont} = strcat(strcat(poly.label,'#V'),num2str(cont));
            end
        end
    else
        if ((nargin == 2) && (iscell(M)))
            if (length(M) > 1)
                error('If the number of vertices and degrees are not specified, the input matrices must not contain the exponent information.');
            end
            M = M{1};
        end
        if (scalar)
            poly.opcode{1} = strcat(poly.label,'#K1');
        elseif ((nargin == 2) && (sum(sum(M)) == 0))
            %Zero matrix
            poly.opcode{1} = strcat(poly.label,'#P1');
        elseif ((nargin == 2) && (size(M,1) == size(M,2)) && (sum(sum(M-eye(length(M)))) == 0))
            %Identity matrix
            poly.opcode{1} = strcat(poly.label,'#P1');
        else
            for cont=1:length(poly.data)
                poly.opcode{cont} = strcat(strcat(poly.label,'#C'),num2str(cont));
            end
        end
    end
    % varargout{1} = poly;
    % if (nargout == 2)
    %     varargout{2} = M;
    % end
    
    poly = rolmipvar(poly);
end


return


function index = gethash(exponent,exptable,jump)
index = 0;
for contsimplex=1:(length(exponent)) 
    if (length(exptable{contsimplex}) > 0)
        [i,j] = find(exptable{contsimplex} - repmat(exponent{contsimplex},size(exptable{contsimplex},1),1) == 0);
        index = index + (find(histc(i,1:size(exptable{contsimplex},1)) == size(exptable{contsimplex},2)) - 1)*jump(contsimplex);
    end
end
index = index + 1;
return