function [resul] = parser_poly(varargin)
% [resul] = parser_poly(expr,vecpoly,newlabel)
%Algorithm to parser the expression to be analyzed and perform
%the calculation.
%vecpoly -> vector with the polynomials
%vecpoly(1).label   -> string with the name of the variable (eg. 'P')
%vecpoly(1).vertices -> number of vertices
%vecpoly(1).data(n) -> data correspondent to the n-th monomial
%vecpoly(1).data(n).exponent -> values of the exponents of the monomial (eg. [0 2])
%vecpoly(1).data(n).value    -> value of the monomial
%
%newlabel: optional, indicates the label received by the result

%Polya relaxation on the polynomial P (example 2 vertices):
%poly_struct(P,'P',vertices,1);
%poly_struct(ones(1,vertices),'1',vertices,1); -->(alpha1 + alpha2)
%[polyrelax] = parser_poly('P*1');

if (nargin == 4)
    expr = varargin{1};
    vecpoly = varargin{2};
    vertices = varargin{3};
    newlabel = varargin{4};
elseif (nargin == 3)
    % The third argument is the number of vertices or the new label
    expr = varargin{1};
    vecpoly = varargin{2};
    if (ischar(varargin{3}))
        newlabel = varargin{3};
        vertices = 0;
    else
        newlabel = [];
        vertices = varargin{3};
    end
elseif (nargin == 2)
    % The second argument is vecpoly or the new label or the vertice
    expr = varargin{1};
    if (ischar(varargin{2}))
        newlabel = varargin{2};
        vecpoly = rolmip('getvar');
        vertices = 0;
    elseif (isnumeric(varargin{2}))
        newlabel = [];
        vecpoly = rolmip('getvar');
        vertices = varargin{2};
    else
        vertices = 0;
        newlabel = [];
        vecpoly = varargin{2};
    end
elseif (nargin == 1)
    expr = varargin{1};
    vecpoly = rolmip('getvar');
    vertices = 0;
    newlabel = [];
end        

%Taking off the space chars
ind = find(expr==' ');
expr(ind) = [];

contexpr = 1;
fila = [];
tamfila = 0;

while (contexpr <= length(expr))
    trans = false; %if the variable is transpose
    labelatual = [];
    subtr = false;
    if (expr(contexpr)=='-')
        %Happens only if the first char is '-' (e.g. -A+B)
        subtr = true;
        contexpr = contexpr + 1;
    end

    %Get the first variable of a product, if it is under parenthesis
    if ((contexpr<=length(expr)) && (expr(contexpr)=='('))
        contexpr = contexpr + 1;
        auxexpr = [];
        stackparenthesis = '(';
        while ((contexpr<=length(expr)) && (~isempty(stackparenthesis)))
            if (expr(contexpr) == ')')
                stackparenthesis(end) = [];
            elseif (expr(contexpr) == '(')
                stackparenthesis = [stackparenthesis '('];
            end
            auxexpr = strcat(auxexpr,expr(contexpr));
            contexpr = contexpr + 1;
        end
        auxexpr(end) = [];
                    
                
        [atual] = parser_poly(auxexpr,vertices);
        vertices = check_vertices(atual,vertices);
        atual.label = strcat('(',atual.label);
        atual.label = strcat(atual.label,')');
        for contcode=1:length(atual.data)
            atual.opcode{contcode} = strcat('(',strcat(atual.opcode{contcode},')'));
        end
        if ((contexpr<=length(expr)) && (expr(contexpr)==char(39))) %transpose
            for conttrans = 1:length(atual.data)
                atual.data(conttrans).value = (atual.data(conttrans).value)';
            end
            atual.label = strcat(atual.label,char(39));
            for contcode=1:length(atual.data)
                atual.opcode{contcode} = strcat(atual.opcode{contcode},char(39));
            end
            contexpr = contexpr + 1;
        end
        if (subtr) %Subtract
            for contsubtr = 1:length(atual.data) 
                atual.data(contsubtr).value = -atual.data(contsubtr).value;
            end
            atual.label = strcat('-',atual.label);
            for contcode=1:length(atual.data)
                atual.opcode{contcode} = strcat('-',atual.opcode{contcode});
            end
        end
    else

        %Get the first variable of a product of variables
        while ((contexpr<=length(expr)) && (expr(contexpr)~='+') && (expr(contexpr)~='-') && (expr(contexpr)~='*'))
            if (expr(contexpr)==char(39)) %transpose
                trans = true;
            else
                labelatual = [labelatual expr(contexpr)];
            end
            contexpr = contexpr + 1;
        end
        
        %Set the number of vertices, if it is not yet set.
        atual = rolmip('getvar',labelatual);
        vertices = check_vertices(atual,vertices);
        
        if (trans) %If it is transpose
            for conttrans=1:length(atual.data)
                atual.data(conttrans).value = (atual.data(conttrans).value)';
            end
            atual.label = strcat(atual.label,char(39));
            for contcode=1:length(atual.data)
                atual.opcode{contcode} = strcat(atual.opcode{contcode},char(39));
            end
        end
        if (subtr) %If the first char is '-'
            for conttrans=1:length(atual.data)
                atual.data(conttrans).value = -atual.data(conttrans).value;
            end
            atual.label = strcat('-',atual.label);
            for contcode=1:length(atual.data)
                atual.opcode{contcode} = strcat('-',atual.opcode{contcode});
            end
        end
        
    end
    %Get the remaining variables of a product of variables, apply the
    %product and store the value into a stack.
    op = '.';
    if (contexpr < length(expr))
        op = expr(contexpr);
        contexpr = contexpr + 1;
        while (op == '*')
            label2 = [];
            trans = false;

            %Parenthesis
            if ((contexpr<=length(expr)) && (expr(contexpr)=='('))
                contexpr = contexpr + 1;
                auxexpr = [];
                stackparenthesis = '(';
                while ((contexpr<=length(expr)) && (~isempty(stackparenthesis)))
                    if (expr(contexpr) == ')')
                        stackparenthesis(end) = [];
                    elseif (expr(contexpr) == '(')
                        stackparenthesis = [stackparenthesis '('];
                    end
                    auxexpr = strcat(auxexpr,expr(contexpr));
                    contexpr = contexpr + 1;
                end
                auxexpr(end) = [];

                [poly2] = parser_poly(auxexpr,vertices);
                vertices = check_vertices(poly2,vertices);
                poly2.label = strcat('(',poly2.label);
                poly2.label = strcat(poly2.label,')');
                for contcode=1:length(poly2.data)
                    poly2.opcode{contcode} = strcat('(',strcat(poly2.opcode{contcode},')'));
                end
                if ((contexpr<=length(expr)) && (expr(contexpr)==char(39))) %transpose
                    for conttrans = 1:length(poly2.data)
                        poly2.data(conttrans).value = (poly2.data(conttrans).value)'; 
                    end
                    poly2.label = strcat(poly2.label,char(39));
                    for contcode=1:length(poly2.data)
                        poly2.opcode{contcode} = strcat(poly2.opcode{contcode},char(39));
                    end
                    contexpr = contexpr + 1;
                end
                
            else
                %No parenthesis
                while ((contexpr<=length(expr)) && (expr(contexpr)~='+') && (expr(contexpr)~='-') && (expr(contexpr)~='*'))
                    if (expr(contexpr)==char(39))
                        trans = true;
                    else
                        label2 = [label2 expr(contexpr)];
                    end
                    contexpr = contexpr + 1;
                end
                
                %Set the number of vertices, if it is not yet set.
                poly2 = rolmip('getvar',label2);
                vertices = check_vertices(poly2,vertices);
                
                if (trans)
                    for conttrans=1:length(poly2.data)
                        poly2.data(conttrans).value = (poly2.data(conttrans).value)';
                    end
                    poly2.label = strcat(poly2.label,char(39));
                    for contcode=1:length(poly2.data)
                        poly2.opcode{contcode} = strcat(poly2.opcode{contcode},char(39));
                    end
                end
                
            end

            [atual] = operation_poly(atual,poly2,op);
            op = '.';
            if (contexpr <= length(expr))
                op = expr(contexpr);
                contexpr = contexpr + 1;
            end
        end
    end
    
    if ((op == '+') || (op == '-'))
        fila{tamfila+2} = op;
    end
    
    fila{tamfila+1} = atual;
    tamfila = tamfila + 2;
end

resul = fila{1};
polyone = rolmip('getvar','1'); %Save the current '1' polynomial
for contsimplex = 1:length(vertices)
    if (vertices(contsimplex) > 0)
        %one = sum_{i=1}^{vertices} alpha_i = 1
        clear M
        degaux = [1 zeros(1,vertices(contsimplex)-1)];
        for cont1 = 1:vertices(contsimplex)
            for cont2 = 1:length(vertices)
                if (cont2 == contsimplex)
                    M{cont1}{cont2} = degaux;
                    degaux = circshift(degaux,[0 1]);
                else
                    M{cont1}{cont2} = [0];
                end
            end
            M{cont1}{cont2+1} = 1;
        end
        
        deg = zeros(1,length(vertices));
        deg(contsimplex) = 1;
        one{contsimplex} = poly_struct(M,'1',vertices,deg);
    end
end
    
if (length(resul.vertices) < length(vertices))
    resul = insert_simplex(resul,vertices(length(resul.vertices)+1:end));
end
contfila = 2;
while (contfila < tamfila)
    op = fila{contfila};
    poly2 = fila{contfila+1};
    %Homogeneize the polynomials to the same degree, if necessary
    if (length(poly2.vertices) < length(vertices))
        poly2 = insert_simplex(poly2,vertices(length(poly2.vertices)+1:end));
    end
    for contsimplex=1:length(vertices)
        deg1 = sum(resul.data(1).exponent{contsimplex});
        deg2 = sum(poly2.data(1).exponent{contsimplex});
        while (deg1 < deg2)
            resul = operation_poly(resul,one{contsimplex},'*');
            deg1 = deg1 + 1;
        end
        while (deg2 < deg1)
            poly2 = operation_poly(poly2,one{contsimplex},'*');
            deg2 = deg2 + 1;
        end
    end
    %Apply the operation
    [resul] = operation_poly(resul,poly2,op);
    contfila = contfila + 2;
end

if (~isempty(polyone)) %Restore the '1' polynomial
    rolmip('setvar',polyone);
end


if (~isempty(newlabel))
    resul.label = newlabel;
end

rolmip('setvar',resul);

return



function vertices = check_vertices(atual,vertices)
%Function to verify and update the current number of vertices

if (sum(vertices) == 0)
    vertices = atual.vertices;
else
    if (length(vertices) > length(atual.vertices))
        ind = find(vertices(1:length(atual.vertices)) == 0);
        vertices(ind) = atual.vertices(ind);
        
        ind = find(vertices(1:length(atual.vertices)).*(atual.vertices) ~= 0);
        if (sum(abs(vertices(ind) - atual.vertices(ind))) ~= 0)
            %Error: different number of vertices!!!!!
            error('The number of vertices of the defined simplexes must always be the same!');
        end
    else
        ind = find(vertices == 0);
        vertices(ind) = atual.vertices(ind);
        if (length(vertices) < length(atual.vertices))
            vertices = [vertices atual.vertices(length(vertices)+1:end)];
        end
        
        ind = find(vertices.*atual.vertices ~= 0);
        if (sum(abs(vertices(ind) - atual.vertices(ind))) ~= 0)
            %Error: different number of vertices!!!!!
            error('The number of vertices of the defined simplexes must always be the same!');
        end
    end
end
    
return