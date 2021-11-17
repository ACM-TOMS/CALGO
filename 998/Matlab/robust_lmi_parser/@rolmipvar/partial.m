function varargout = partial(varargin)
% Procedure to calculate the partial derivate over a simplex of a
% parameter-dependent matrix.
%
% [dpoly] = partial(poly,simpnum) calculates the partial derivative of the
% polynomial poly over a given simplex. The input simpnum describes the
% number of the simplex over which the derivative must be performed, and
% the output dpoly is a cell structure of N indexes, being N the number of
% vertices of the simplex described by simpnum.
%
% [dpoly] = partial(poly) calculates the partial derivative over the first
% simplex.
%
% [dpoly] = partial(poly,labelout,simpnum) also defines the label of
% the output polynomial as the string labelout. If not informed, the 
% output label, as standard, adds the string '_dX#' to the input label,
% being X the letter used to describe the simplex and # the number of the 
% vertex.
%
% [dpoly] = partial(poly,labelout) defines the label of the output
% polynomial as the string labelout. The partial derivative is calculated
% over the first simplex.
%
% Author: Cristiano M. Agulhari
% 2017, Nov, 10

if ((nargin < 1) || (nargin > 3))
    error('Input error. Type ''help partial'' for more details');
    return
else
    polyin = varargin{1};
    if (nargin == 1)
        simpnum = 1;
        labelout = strcat('part_',strcat(char('a'+simpnum-1),strcat('_',polyin.label)));
    elseif (nargin == 2)
        aux = varargin{2};
        if (isstr(aux)) %String
            simpnum = 1;
            labelout = aux;
        else
            simpnum = varargin{2};
            labelout = strcat('part_',strcat(char('a'+simpnum-1),strcat('_',polyin.label)));
        end
    elseif (nargin == 3)
        labelout = varargin{2};
        simpnum = varargin{3};
    end
end

if (sum(polyin.vertices) == 0)
    %In this case, the polynomial is a constant matrix
    varargout{1} = rolmipvar(zeros(size(polyin.data(1).value)),labelout);
    return
end

numsimplexes = length(polyin.vertices);
if (simpnum > numsimplexes)
    error('Index exceeds the number of simplexes.');
    return
end

for contsimplex=1:numsimplexes
    degree(contsimplex) = sum(polyin.data(1).exponent{contsimplex});
end

if (sum(degree) == 0)
    %In this case, the polynomial is a constant matrix
    varargout{1} = rolmipvar(zeros(size(polyin.data(1).value)),labelout);
    return
end

if (degree(simpnum) == 0)
    %The polynomial does not depend on the simplex variable
    varargout{1} = rolmipvar(zeros(size(polyin.data(1).value)),labelout);
    return
end

% Define the vertices and degrees of the simplexes related to the
% derivatives
vertices = polyin.vertices;

degree(simpnum) = degree(simpnum)-1;
[jump,exptable] = create_hash_table(vertices,degree);

% Compose the polynomial structure of the derivative
for contvert=1:vertices(simpnum)
    poly{contvert}.label = labelout;
    poly{contvert}.vertices = vertices;
    poly{contvert}.data = [];
    poly{contvert}.opcode = [];
    poly{contvert}.opcodein = [];
end

for cont=1:length(polyin.data)
    exporig = polyin.data(cont).exponent;
    if (sum(exporig{simpnum}) > 0)
        expnow = exporig{simpnum};
        for contexp=1:length(expnow)
            if (expnow(contexp) > 0)
                exponents = exporig;
                val = exponents{simpnum}(contexp);
                exponents{simpnum}(contexp) = exponents{simpnum}(contexp) - 1;
                indresul = gethash(exponents,exptable,jump);
                if ((indresul > length(poly{contexp}.data)) || (~isfield(poly{contexp}.data(indresul),'value')) || (isempty(poly{contexp}.data(indresul).value)))
                    poly{contexp}.data(indresul).value = 0;
                end
                poly{contexp}.data(indresul).exponent = exponents;
                poly{contexp}.data(indresul).value = poly{contexp}.data(indresul).value + (exponents{simpnum}(contexp)+1)*polyin.data(cont).value;
                
                %Update the opcode and opcodein
                if ~((indresul > length(poly{contexp}.opcodein)) || (isempty(poly{contexp}.opcodein{indresul})))
                    %If the opcodein already exists
                    poly{contexp}.opcodein{indresul} = strcat(poly{contexp}.opcodein{indresul},'+');
                else
                    poly{contexp}.opcodein{indresul} = [];
                end
                if ((indresul > length(poly{contexp}.opcode)) || (isempty(poly{contexp}.opcode{indresul})))
                    poly{contexp}.opcode{indresul} = strcat(labelout,strcat('#T',num2str(indresul)));
                end
                if (val~=1)
                    poly{contexp}.opcodein{indresul} = strcat(strcat(poly{contexp}.opcodein{indresul},num2str(val)),'*');
                end
%                 poly.opcodein{indresul} = strcat(poly.opcodein{indresul},'H[');
%                 poly.opcodein{indresul} = strcat(strcat(strcat(strcat(poly.opcodein{indresul},strcat(strcat(num2str(contexp),','),num2str(contderiv))),']#H'),num2str(contsimplex)),'*');
                poly{contexp}.opcodein{indresul} = strcat(poly{contexp}.opcodein{indresul},polyin.opcode{cont});
            end
        end
    end
end
                    
for contvert=1:vertices(simpnum)
    output{contvert} = rolmipvar(poly{contvert});
end
varargout{1} = output;

return


%function [indhash,base] = create_hash_table(vertices,degree)
function [jump,exptable] = create_hash_table(vertices,degree)
jump = 1;
numelem = 1;
for contsimplex=1:length(vertices)
    [exptable{contsimplex}] = generate_homogeneous_exponents(vertices(contsimplex),degree(contsimplex));
    if (contsimplex > 1)
        jump(contsimplex) = numelem;
    end
    numelem = numelem*size(exptable{contsimplex},1);
end
%         
% numelem = 1;
% vertnow = 0;
% base = [];
% for contsimplex=1:length(vertices)
%     [exponents{contsimplex}] = generate_homogeneous_exponents(vertices(contsimplex),degree(contsimplex));
%     
%     if (degree(contsimplex) > 0)
%         if (sum(degree) > 1)
%             base = [base (sum(degree)*ones(1,vertices(contsimplex))).^(vertnow:vertnow+vertices(contsimplex)-1)];
%             vertnow = vertnow + vertices(contsimplex);
%         else
%             base = [base (2*ones(1,vertices(contsimplex))).^(vertnow:vertnow+vertices(contsimplex)-1)];
%             vertnow = vertnow + vertices(contsimplex);
%         end
%     end
%     numelem = numelem*size(exponents{contsimplex},1);
% end
% 
% contexp = ones(1,length(vertices)+1);
% for conttotal=1:numelem
%     vetexponent = [];
%     for contsimplex=1:length(exponents)
%         polyexponent{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);
%         if (sum(polyexponent{contsimplex}) > 0)
%             vetexponent = [vetexponent polyexponent{contsimplex}];
%         end
%     end
%     indhash(conttotal) = sum(base.*vetexponent);
%     
%     contexp(1) = contexp(1) + 1;
%     aux = 1;
%     while ((aux <= length(vertices)) && (contexp(aux) > size(exponents{aux},1)))
%         contexp(aux) = 1;
%         contexp(aux+1) = contexp(aux+1) + 1;
%         aux = aux + 1;
%     end
% end
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
% function index = gethash(exponent,hash,base)
% vetexponent = [];
% for contsimplex=1:(length(exponent)) 
%    if (sum(exponent{contsimplex}) > 0)
%         vetexponent = [vetexponent exponent{contsimplex}];
%    end
% end
% 
% index = find(sum(base.*(vetexponent))==hash);
% return
