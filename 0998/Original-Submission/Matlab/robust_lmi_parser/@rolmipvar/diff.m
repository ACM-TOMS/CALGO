function varargout = diff(varargin)
% Procedure to calculate the time-derivate of a time-varying
% parameter-dependent matrix.
%
% [dotpoly] = diff(poly,dotbounds) calculates the
% derivative of the polynomial poly, being the bounds of the
% variation rates of the parameters given by the rows of the vector
% dotbounds. For example, considering that the derivative of a1 is bounded
% by x <= dot(a1) <= y and of a2 is bounded by w <= dot(p2) <= z, then
% dotbounds = [x y; w z]. 
% If the variable is a multi-simplex, then dotbounds{i} is the vector with
% the variation rates bounds of the parameters in the i-th simplex.
%
% [dotpoly] = diff(poly,labelout,dotbounds) also defines the label of
% the output polynomial as the string labelout. If not informed, the 
% output label, as standard, adds the string 'dot_' to the input label.
%
% [dotpoly] = diff(poly,bounds,dotbounds) is used when the
% polynomial stems from an affine representation in the form  m(p) = m0 +
% p1*m1 + ... + pn*mn, being the bounds of the parameters given by the
% vector bounds. In this case, the vector dotbounds contains the bounds of
% the derivatives of the parameters p.
%
% [dotpoly] = diff(poly,labelout,bounds,dotbounds) also defines the label 
% of the output polynomial as the string labelout.

% Author: Cristiano M. Agulhari
% 2016, Feb, 3

if ((nargin < 2) || (nargin > 4))
    error('Input error. Type ''help diff'' for more details');
    return
else
    polyin = varargin{1};
    if (nargin == 2)
        dotbounds = varargin{2};
        labelout = strcat('dot_',polyin.label);
    elseif (nargin == 3)
        if (isstr(varargin{2}))
            dotbounds = varargin{3};
            labelout = varargin{2};
        else
            bounds = varargin{2};
            dotbounds = varargin{3};
            dotbounds = affine2simplex(bounds,dotbounds);
            labelout = strcat('dot_',polyin.label);
        end
    elseif (nargin == 4)
        bounds = varargin{3};
        dotbounds = varargin{4};
        dotbounds = affine2simplex(bounds,dotbounds);
        labelout = varargin{2};
    end
end

if (~iscell(dotbounds))
    aux = dotbounds;
    clear dotbounds;
    dotbounds{1} = aux;
    clear aux;
end

if (sum(polyin.vertices) == 0)
    %In this case, the polynomial is a constant matrix
    varargout{1} = rolmipvar(zeros(size(polyin.data(1).value)),labelout);
    return
end

numsimplexes = length(polyin.vertices);
for contsimplex=1:numsimplexes
    degree(contsimplex) = sum(polyin.data(1).exponent{contsimplex});
end

if (sum(degree) == 0)
    %In this case, the polynomial is a constant matrix
    varargout{1} = rolmipvar(zeros(size(polyin.data(1).value)),labelout);
    return
end

% Define the vertices and degrees of the simplexes related to the
% derivatives
vertices = polyin.vertices;
for contsimplex=1:numsimplexes
    H{contsimplex} = matrixH(dotbounds{contsimplex});
    vertices = [vertices size(H{contsimplex},2)];
    degree = [degree 1];
end  

[jump,exptable] = create_hash_table(vertices,degree);

% Compose the polynomial structure of the derivative
poly.label = labelout;
poly.vertices = vertices;
poly.data = [];
poly.opcode = [];
poly.opcodein = [];
for cont=1:length(polyin.data)
    exporig = polyin.data(cont).exponent;
    for aux=numsimplexes+1:2*numsimplexes
        %exporig{aux} = [1 zeros(1,vertices(aux)-1)];
        exporig{aux} = zeros(1,vertices(aux));
    end
    for contsimplex=1:numsimplexes
        expnow = exporig{contsimplex};
        for contexp=1:length(expnow)
            %Derivate the simplex variable related to contexp
            if (expnow(contexp) > 0)
                exponents = exporig;
                val = exponents{contsimplex}(contexp);
                
                exponents{contsimplex}(contexp) = exponents{contsimplex}(contexp) - 1;
                %Applying the derivative on the simplex variable related to
                %exponents{contsimplex}(contexp)
                
                %First step: Homogeneize the current simplex variable
                for contatual=1:length(exponents{contsimplex})
                    exponents{contsimplex}(contatual) = exponents{contsimplex}(contatual) + 1;
                    indderiv = contsimplex + numsimplexes; %index of the derivative simplex
                    for contderiv=1:length(exponents{indderiv})
                        exponents{indderiv}(contderiv) = 1;
                        %Second step: Homogeneize on all the derivative
                        %simplexes, except the current derivative
                        indderivhomog = numsimplexes+1:2*numsimplexes;
                        indderivhomog(contsimplex) = [];
                        
                        for contderivhomog=1:length(indderivhomog)
                            %Prepare the exponents
                            exponents{indderivhomog(contderivhomog)}(1) = 1;
                        end
                        finish = false;
                        while (~finish)
                            indresul = gethash(exponents,exptable,jump);
                            if ((indresul > length(poly.data)) || (~isfield(poly.data(indresul),'value')) || (isempty(poly.data(indresul).value)))
                                poly.data(indresul).value = 0;
                            end
                            poly.data(indresul).value = poly.data(indresul).value + val*H{contsimplex}(contexp,contderiv)*polyin.data(cont).value;
                            poly.data(indresul).exponent = exponents;
                            
                             %Update the opcode and opcodein
                            if ~((indresul > length(poly.opcodein)) || (isempty(poly.opcodein{indresul})))
                                %If the opcodein already exists
                                poly.opcodein{indresul} = strcat(poly.opcodein{indresul},'+');
                            else
                                poly.opcodein{indresul} = [];
                            end
                            if ((indresul > length(poly.opcode)) || (isempty(poly.opcode{indresul})))
                                poly.opcode{indresul} = strcat(labelout,strcat('#D',num2str(indresul)));
                            end
                            if (val~=1)
                                poly.opcodein{indresul} = strcat(strcat(poly.opcodein{indresul},num2str(val)),'*');
                            end
                            poly.opcodein{indresul} = strcat(poly.opcodein{indresul},'H[');
                            poly.opcodein{indresul} = strcat(strcat(strcat(strcat(poly.opcodein{indresul},strcat(strcat(num2str(contexp),','),num2str(contderiv))),']#H'),num2str(contsimplex)),'*');
                            poly.opcodein{indresul} = strcat(poly.opcodein{indresul},polyin.opcode{cont});
                            
                            isok = false;
                            if (isempty(indderivhomog))
                                isok = true;
                                finish = true;
                            end
                            contderivhomog = 1;
                            while (~isok)
                                exponents{indderivhomog(contderivhomog)} = circshift(exponents{indderivhomog(contderivhomog)},[0 1]);
                                if (exponents{indderivhomog(contderivhomog)}(1) == 1)
                                    contderivhomog = contderivhomog + 1;
                                    if (contderivhomog > length(indderivhomog))
                                        isok = true;
                                        finish = true; %Finished the homogenization
                                    end
                                else
                                    isok = true;
                                end
                            end
                        end 
                        exponents{indderiv}(contderiv) = 0;
                    end
                    exponents{contsimplex}(contatual) = exponents{contsimplex}(contatual) - 1;
                end
            end
        end
    end
end
                    
varargout{1} = rolmipvar(poly);

return

function dotbsimplex = affine2simplex(bounds,dotbounds)
% Function that receives the bounds and dotbounds of an affine
% representation and returns the bounds of the variation rates in the
% simplex representation (dotbsimplex).

for cont=1:size(dotbounds,1)
    auxmin = dotbounds(cont,1)/(bounds(cont,1) - bounds(cont,2));
    auxmax = dotbounds(cont,2)/(bounds(cont,1) - bounds(cont,2));
    if (auxmin < auxmax)
        dotbsimplex{cont} = [auxmin auxmax; -auxmax -auxmin];
    else
        dotbsimplex{cont} = [auxmax auxmin; -auxmin -auxmax];
    end
end
return


function H = matrixH(dotbounds)
%Procedure to determine the H matrix for the derivative of polynomial
%matrices

if (exist('Polyhedron')) %If MPT Toolbox is installed
    N = size(dotbounds,1);
    A = kron(eye(N),[1;-1]);
    M = [dotbounds(:,2) -dotbounds(:,1)];
    b = reshape(M',size(dotbounds,1)*size(dotbounds,2),1);
    Ae = ones(1,N);
    be = 0;
    
    
    P = Polyhedron('A',A,'b',b,'Ae',Ae,'be',be);
    H = P.V';
else
    
    
    N = size(dotbounds,1);
    H = [];
    for i=1:N
        Y{i}=dotbounds(i,:)';
        Z{i}=dotbounds';
        Z{i}(:,i)=[];
        
        
        C=[];
        for j=1:N-1
            C=[C; kron(ones(1,2^(j-1)),kron(Z{i}(:,j)',ones(1,2^(N-1-j))))];
        end
        if size(C,1)==1
            V=-C;
        else
            V=-sum(C);
        end
        for j=1:2^(N-1)
            if V(j) >= Y{i}(1) && V(j) <= Y{i}(2)%if V(j) >= Y{i}(1) & V(j) <= Y{i}(2)
                c=C(:,j);
                col = [c(1:i-1,1); V(j); c(i:end,1)];
                
                if i ==1;
                    H=[H col];
                end
                count=0;
                for r=1:size(H,2);
                    v=H(:,r);
                    if roundn(v-col,-8) == zeros(N,1);
                        count=count+1;
                    end
                end
                
                if count ==0;
                    H=[H col];
                end
                
            end
        end
    end
end

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