function varargout = lmifiles(varargin)
% Create and manipulate a .m file to solve faster a given set of LMIs.
%
% fid = lmifiles('open',filename);
% fid = lmifiles('o',filename);
% Open the .m file with name given by 'filename'
%
% lmifiles('insert',fid,Term,ineq);
% lmifiles('i',fid,Term,ineq);
% Insert a set of LMIs in the .m file with descriptor 'fid'. The LMIs are
% defined in the structure 'Term' and the signal of the inequality is 
% given by 'ineq'.
%
% lmifiles('i', fid, label_obj)
% Insert the variable whose label is label_obj as the minimizing
% objective function.
%
% lmifiles('i',fid,'sdpsettings','NAME1',VALUE1,'NAME2',VALUE2,...)
% Define the settings used in the solvesdp command to solve the LMIs.
%
% lmifiles('close',fid);
% lmifiles('c',fid);
% Close the .m file with descriptor 'fid'.


action = varargin(1);

if (strcmp(action,'open') || strcmp(action,'o'))
    %Open the file
    filename = varargin(2);
    filenametemp = strcat(filename,'.temp');
    fid = fopen(char(filenametemp),'wt');
    fprintf(fid,'LMIs = [];\n\n');
    variables{1024} = [];
    numLMI = 1;
    indexes = [];
    save lmifilestemp variables numLMI indexes filename %Used to track the variables and constants
    varargout{1} = fid;
    fidaux = fopen('varaux_temporary_123.temp','wt');
    rolmip('setauxfile',fidaux);
    rolmip('setobj',[]);
end


if (strcmp(action,'insert') || strcmp(action,'i'))
    %Insert a LMI term on the file
    if (strcmp(varargin{3},'sdpsettings'))
        fid = varargin{2};
        opt = [];
        for cont=4:length(varargin)
            if (cont > 4)
                opt = strcat(opt,',');
            end

            if (isstr(varargin{cont}))
                opt = strcat(opt,'''');
                opt = strcat(opt,varargin{cont});
                opt = strcat(opt,'''');
            else
                opt = strcat(opt,num2str(varargin{cont}));
            end
        end
        rolmip('setopt',opt);
        Term = [];
    elseif (nargin == 6)
        fid = varargin{2};
        Term = varargin{3};
        vecpoly = varargin{4};
        vertices = varargin{5};
        ineq = varargin{6};
    elseif (nargin == 5)
        fid = varargin{2};
        Term = varargin{3};
        vecpoly = varargin{4};
        ineq = varargin{5};
        vertices = 0;
        conti = 1;
        contj = 1;
        while ((sum(vertices) == 0) && (conti <= size(Term,1)) && (contj <= size(Term,2)))
            if (~isempty(Term{conti,contj}))
                vertices = Term{conti,contj}.vertices;
            else
                vertices = 0;
            end
            contj = contj + 1;
            if (contj > size(Term,2))
                conti = conti + 1;
                contj = 1;
            end
        end
    elseif (nargin == 4)
        fid = varargin{2};
        Term = varargin{3};
        ineq = varargin{4};
        vecpoly = rolmip('getvar');
        vertices = 0;
        conti = 1;
        contj = 1;
        while ((sum(vertices) == 0) && (conti <= size(Term,1)) && (contj <= size(Term,2)))
            if (~isempty(Term{conti,contj}))
                vertices = Term{conti,contj}.vertices;
            else
                vertices = 0;
            end
            contj = contj + 1;
            if (contj > size(Term,2))
                conti = conti + 1;
                contj = 1;
            end
        end
    elseif (nargin == 3)
        fid = varargin{2};
        obj = varargin{3};
        rolmip('setobj',obj);
        Term = [];
    end

    if (~isempty(Term))
        Term = construct_lmi_terms(Term,vertices,'Term');

        update_variables(Term,vecpoly); %Update the variables and constants found in Term

        write_matrix_file(fid,Term,ineq); %Writes 'Term' in the file and sets the LMI
        
    end
end


if (strcmp(action,'close') || strcmp(action,'c'))
    load lmifilestemp
    
    fidtemp = varargin{2};
    fclose(fidtemp);
    fidtemp = fopen(char(strcat(filename,'.temp')),'r');
    
    fid = fopen(char(strcat(filename,'.m')),'wt');
    
    %Function name and input arguments
    fprintf(fid,'function output = %s(',char(filename));
    contC = 0; contV = 0; contP = 0; contK = 0; contS = 0; contA = 0; 
    contZ = 0; contF = 0; contD = 0;
    %Type C: constant polynomial
    %Type V: variable polynomial
    %Type P: predefinded identity or zero matrices
    %Type K: scalar constant
    %Type S: scalar variable
    %Type A: auxiliar matrix (created from construct_lmi)
    %Type Z: affine constant matrix
    %Type F: affine variable matrix
    %Type H: matrix H, generated from the bounds of derivatives of params.
    %Type D: polynomials dependent on other polynomials
    for cont=1:length(variables)
        if (~isempty(variables{cont}))
            if (variables{cont}.type == 'C')
                contC = contC + 1;
                varC{contC} = variables{cont};
                if (contK + contC + contZ > 1)
                    fprintf(fid,',');
                end
                fprintf(fid,'%s',variables{cont}.name);
            elseif (variables{cont}.type == 'V')
                contV = contV + 1;
                varV{contV} = variables{cont};
            elseif (variables{cont}.type == 'P')
                contP = contP + 1;
                varP{contP} = variables{cont};
            elseif (variables{cont}.type == 'K')
                contK = contK + 1;
                varK{contK} = variables{cont};
                if (contK + contC + contZ > 1)
                    fprintf(fid,',');
                end
                fprintf(fid,'%s',variables{cont}.name);
            elseif (variables{cont}.type == 'S')
                contS = contS + 1;
                varS{contS} = variables{cont};
            elseif (variables{cont}.type == 'A')
                contA = contA + 1;
                varA{contA} = variables{cont};
            elseif (variables{cont}.type == 'Z')
                contZ = contZ + 1;
                varZ{contZ} = variables{cont};
                if (contK + contC + contZ > 1)
                    fprintf(fid,',');
                end
                fprintf(fid,'%s',variables{cont}.name);
            elseif (variables{cont}.type == 'F')
                contF = contF + 1;
                varF{contF} = variables{cont};
            elseif (variables{cont}.type == 'D')
                contD = contD + 1;
                varD{contD} = variables{cont};
            end
        end
    end
    if (contZ > 0)
        fprintf(fid,',bounds)\n');
        fprintf(fid,'%% bounds: The bounds of the parameters of the affine representations.\n\n');
    else
        fprintf(fid,')\n\n');
    end
    

    %Obtain the dimensions of the constants #C
    if (contC > 0)
        for cont=1:length(varC)
            rowvar = strcat('row',varC{cont}.name);
            colvar = strcat('col',varC{cont}.name);
            %Detects if it is a cell array or a concatenation
            fprintf(fid,'if (iscell(%s))\n',varC{cont}.name);
            fprintf(fid,'if (iscell(%s{1}))\n',varC{cont}.name);
            fprintf(fid,'[%s, %s] = size(%s{1}{length(%s{1})});\n',rowvar,colvar,varC{cont}.name,varC{cont}.name);
            fprintf(fid,'%s = getcoefsmonomials(%s);\n',varC{cont}.name,varC{cont}.name);
            fprintf(fid,'else\n');
            fprintf(fid,'[%s, %s] = size(%s{1});\n',rowvar,colvar,varC{cont}.name);
            fprintf(fid,'end\n');
            fprintf(fid,'else\n');
            fprintf(fid,'[%s, %s] = size(%s);\n',rowvar,colvar,varC{cont}.name);
            fprintf(fid,'%s = %s/%d;\n',colvar,colvar,varC{cont}.nummon);
            fprintf(fid,'for cont = 1:%d\n',varC{cont}.nummon);
            fprintf(fid,'%s{cont} = %s(:,(cont-1)*%s+1:cont*%s);\n',strcat(varC{cont}.name,'lmiaux'),varC{cont}.name,colvar,colvar);
            fprintf(fid,'end\n'); 
            fprintf(fid,'%s = %s;\n',varC{cont}.name,strcat(varC{cont}.name,'lmiaux'));
            fprintf(fid,'clear %s;\n',strcat(varC{cont}.name,'lmiaux'));
            fprintf(fid,'end\n\n');
        end
    end

    fprintf(fid,'\n');

    
    %Obtain the dimensions of the affine constants #Z
    if (contZ > 0)
        for cont=1:length(varZ)
            rowvar = strcat('row',varZ{cont}.name);
            colvar = strcat('col',varZ{cont}.name);
            %Transform from the affine representation to the simplex
            fprintf(fid,'if (iscell(%s))\n',varZ{cont}.name);
            fprintf(fid,'[%s, %s] = size(%s{1});\n',rowvar,colvar,varZ{cont}.name);
            fprintf(fid,'else\n');
            fprintf(fid,'[%s, %s] = size(%s);\n',rowvar,colvar,varZ{cont}.name);
            fprintf(fid,'%s = %s/%d;\n',colvar,colvar,varZ{cont}.nummon);
            fprintf(fid,'for cont = 1:%d\n',varZ{cont}.nummon);
            fprintf(fid,'%s{cont} = %s(:,(cont-1)*%s+1:cont*%s);\n',strcat(varZ{cont}.name,'lmiaux'),varZ{cont}.name,colvar,colvar);
            fprintf(fid,'end\n'); 
            fprintf(fid,'%s = %s;\n',varZ{cont}.name,strcat(varZ{cont}.name,'lmiaux'));
            fprintf(fid,'clear %s;\n',strcat(varZ{cont}.name,'lmiaux'));
            fprintf(fid,'end\n');
            fprintf(fid,'%s = convertaffine(%s,bounds);\n\n',varZ{cont}.name,varZ{cont}.name);
        end
    end

    fprintf(fid,'\n');
   
    %Define variables type V
    if (contV > 0)
        for cont=1:length(varV)
            if (varV{cont}.nummon == 1)
                fprintf(fid,'%s{1} = sdpvar(%s,%s,''%s'');\n',varV{cont}.name,varV{cont}.dimension{1},varV{cont}.dimension{2},varV{cont}.varstruct);
            else
                fprintf(fid,'for cont=1:%d\n',varV{cont}.nummon);
                fprintf(fid,'%s{cont} = sdpvar(%s,%s,''%s'');\n',varV{cont}.name,varV{cont}.dimension{1},varV{cont}.dimension{2},varV{cont}.varstruct);
                fprintf(fid,'end\n');
            end
        end
    end

    fprintf(fid,'\n');
    
    %Define variables type F
    if (contF > 0)
        for cont=1:length(varF)
            if (varF{cont}.nummon == 1)
                fprintf(fid,'%s{1} = sdpvar(%s,%s,''%s'');\n',varF{cont}.name,varF{cont}.dimension{1},varF{cont}.dimension{2},varF{cont}.varstruct);
            else
                fprintf(fid,'for cont=1:%d\n',varF{cont}.nummon);
                fprintf(fid,'%s{cont} = sdpvar(%s,%s,''%s'');\n',varF{cont}.name,varF{cont}.dimension{1},varF{cont}.dimension{2},varF{cont}.varstruct);
                fprintf(fid,'end\n');
                fprintf(fid,'%s = convertaffine(%s,bounds);\n\n',varF{cont}.name,varF{cont}.name);
            end
        end
    end

    fprintf(fid,'\n');

    %Define variable scalars type S
    if (contS > 0)
        for cont=1:length(varS)
            fprintf(fid,'%s = sdpvar(1,1);\n',varS{cont}.name);
        end
    end

    fprintf(fid,'\n');

    %Set predefined variables P
    if (contP > 0)
        for cont=1:length(varP)
            if (strcmp(varP{cont}.varstruct,'identity'))
                fprintf(fid,'%s = eye(%s);\n',varP{cont}.name,varP{cont}.dimension{1});
            else
                fprintf(fid,'%s = zeros(%s,%s);\n',varP{cont}.name,varP{cont}.dimension{1},varP{cont}.dimension{2});
            end
        end
    end

    fprintf(fid,'\n');
    
    %Set the auxiliar variables A
    fidaux = [];
    if (contA > 0)
        %Print the content of the auxiliar temporary file
        fidaux = rolmip('getauxfile');
        fclose(fidaux);
        fidaux = fopen('varaux_temporary_123.temp','r');
        rolmip('setauxfile',fid);
        tline = fgets(fidaux);
        while (tline ~= -1)
            fprintf(fid,tline);
            tline = fgets(fidaux);
        end
    end
    
    fprintf(fid,'\n\n');
    
    %Set the dependent variables D
    if (contD > 0)
        for cont=1:length(varD)
            
            polaux = rolmip('getvar',varD{cont}.name);
            fprintf(fid,'\n');
            for contvar=1:length(polaux.opcode)
                expr = polaux.opcodein{contvar};
                toprint = strcat(varD{cont}.name,strcat('{',strcat(num2str(contvar),'} = ')));
                contexpr = 1;
                while (contexpr <= length(expr))
                    varnow = [];
                    while (~strcmp(expr(contexpr),'#'))
                        varnow = strcat(varnow,expr(contexpr));
                        contexpr = contexpr + 1;
                    end
                    contexpr = contexpr + 2;
                    varnow = strcat(varnow,'{');
                    while ((contexpr <= length(expr)) && ~strcmp(expr(contexpr),'+') && ~strcmp(expr(contexpr),'-') && ~strcmp(expr(contexpr),'*') && ~strcmp(expr(contexpr),char(39)) && ~strcmp(expr(contexpr),')') && ~strcmp(expr(contexpr),'('))
                        varnow = strcat(varnow,expr(contexpr));
                        contexpr = contexpr + 1;
                    end
                    varnow = strcat(varnow,'}');
                    while ((contexpr <= length(expr)) && (strcmp(expr(contexpr),'+') || strcmp(expr(contexpr),'-') || strcmp(expr(contexpr),'*') || strcmp(expr(contexpr),char(39)) || strcmp(expr(contexpr),')') || strcmp(expr(contexpr),'(')))
                        varnow = strcat(varnow,expr(contexpr));
                        contexpr = contexpr + 1;
                    end
                    toprint = strcat(toprint,varnow);
                end
                toprint = strcat(toprint,'; \n');
                fprintf(fid,toprint);
            end
        end
    end
        
    

    %Insert the LMIs in the file
    while (~feof(fidtemp))
        rowtemp = fgets(fidtemp);
        if (double(rowtemp(end)) == 10)
            rowtemp(end) = [];
        end
        fprintf(fid,'%s',rowtemp);
    end
    
    fprintf(fid,'\n\n');
    
    obj = rolmip('getobj');
    if (isempty(obj))
        obj = '[]';
    end
    
    opt = rolmip('getopt');
    
    %fprintf(fid,'solvesdp(LMIs,%s,sdpsettings(''solver'',''sedumi'',''verbose'',0));\n',obj);
    fprintf(fid,'solvesdp(LMIs,%s,sdpsettings(%s));\n',obj,opt);
    
    fprintf(fid,'[primal1,dual1] = checkset(LMIs);\n');
    
    fprintf(fid,'if (sum(primal1+1e-7 < 0) == 0)\n');
    fprintf(fid,'output.feas = 1;\n\n');
    if (contV > 0)
        for cont=1:length(varV)
            if (varV{cont}.nummon == 1)
                fprintf(fid,'output.%s{1} = double(%s{1});\n',varV{cont}.name,varV{cont}.name);
            else
                fprintf(fid,'for cont=1:%d\n',varV{cont}.nummon);
                fprintf(fid,'output.%s{cont} = double(%s{cont});\n',varV{cont}.name,varV{cont}.name);
                fprintf(fid,'end\n');
            end
        end
    end
    
    if (contF > 0)
        for cont=1:length(varF)
            if (varF{cont}.nummon == 1)
                fprintf(fid,'output.%s{1} = double(%s{1});\n',varF{cont}.name,varF{cont}.name);
            else
                fprintf(fid,'for cont=1:length(%s)\n',varF{cont}.name);
                fprintf(fid,'output.%s{cont} = double(%s{cont});\n',varF{cont}.name,varF{cont}.name);
                fprintf(fid,'end\n');
            end
        end
    end
    
    if (contS > 0)
        for cont=1:length(varS)
            fprintf(fid,'output.%s = double(%s);\n',varS{cont}.name,varS{cont}.name);
        end
    end
    
    fprintf(fid,'else\n');
    fprintf(fid,'output.feas = 0;\n');
    fprintf(fid,'end\n\n');
            
    fprintf(fid,'return');
    
    %Now print the function getcoefsmonomials, responsible for converting
    %the representation A{1} = {[1 0],A1} to A{1} = A1.
    
    fprintf(fid,'\n\n\n function R = getcoefsmonomials(M)\n');

    fprintf(fid,'numvertices = length(M{1})-1;\n');
    fprintf(fid,'for cont = 1:numvertices\n');
    fprintf(fid,'vertices(cont) = length(M{1}{cont});\n');
    fprintf(fid,'degree(cont) = sum(M{1}{cont});\n');
    fprintf(fid,'end\n');

    fprintf(fid,'numelem = 1;\n');
    fprintf(fid,'vertnow = 0;\n');
    fprintf(fid,'base = [];\n');
    fprintf(fid,'for contsimplex=1:length(vertices)\n');
    fprintf(fid,'[exponents{contsimplex}] = generate_homogeneous_exponents(vertices(contsimplex),degree(contsimplex));\n');

    fprintf(fid,'if (degree(contsimplex) > 0)\n');
    fprintf(fid,'if (sum(degree) > 1)\n');
    fprintf(fid,'base = [base (sum(degree)*ones(1,vertices(contsimplex))).^(vertnow:vertnow+vertices(contsimplex)-1)];\n');
    fprintf(fid,'vertnow = vertnow + vertices(contsimplex);\n');
    fprintf(fid,'else\n');
    fprintf(fid,'base = [base (2*ones(1,vertices(contsimplex))).^(vertnow:vertnow+vertices(contsimplex)-1)];\n');
    fprintf(fid,'vertnow = vertnow + vertices(contsimplex);\n');
    fprintf(fid,'end\n');
    fprintf(fid,'end\n');
    fprintf(fid,'numelem = numelem*size(exponents{contsimplex},1);\n');
    fprintf(fid,'end\n');

    fprintf(fid,'contexp = ones(1,length(vertices)+1);\n');
    fprintf(fid,'for conttotal=1:numelem\n');
    fprintf(fid,'vetexponent = [];\n');
    fprintf(fid,'for contsimplex=1:length(exponents)\n');
    fprintf(fid,'polyexponent{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);\n');
    fprintf(fid,'if (sum(polyexponent{contsimplex}) > 0)\n');
    fprintf(fid,'vetexponent = [vetexponent polyexponent{contsimplex}];\n');
    fprintf(fid,'end\n');
    fprintf(fid,'end\n');
    fprintf(fid,'indhash(conttotal) = sum(base.*vetexponent);\n');

    fprintf(fid,'contexp(1) = contexp(1) + 1;\n');
    fprintf(fid,'aux = 1;\n');
    fprintf(fid,'while ((aux <= length(vertices)) && (contexp(aux) > size(exponents{aux},1)))\n');
    fprintf(fid,'contexp(aux) = 1;\n');
    fprintf(fid,'contexp(aux+1) = contexp(aux+1) + 1;\n');
    fprintf(fid,'aux = aux + 1;\n');
    fprintf(fid,'end\n');
    fprintf(fid,'end\n');

    fprintf(fid,'for cont=1:length(M)\n');
    fprintf(fid,'for contsimplex=1:(length(M{cont})-1)\n');
    fprintf(fid,'exponent{contsimplex} = M{cont}{contsimplex};\n');
    fprintf(fid,'end\n');

    fprintf(fid,'indresul = gethash(exponent,indhash,base);\n');
    fprintf(fid,'R{indresul} = M{cont}{length(M{cont})};\n');
    fprintf(fid,'end\n');

    fprintf(fid,'%Insert the zero-monomials\n');

    fprintf(fid,'contexp = ones(1,length(vertices)+1);\n');
    fprintf(fid,'for conttotal=1:numelem\n');
    fprintf(fid,'exponent = [];\n');
    fprintf(fid,'for contsimplex=1:(length(M{cont})-1)\n');
    fprintf(fid,'exponent{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);\n');
    fprintf(fid,'end\n');

    fprintf(fid,'indresul = gethash(exponent,indhash,base);\n');
    fprintf(fid,'if ((indresul > length(R)) || (isempty(R{indresul})))\n');
    fprintf(fid,'R{indresul} = zeros(size(M{1}{length(M{1})}));\n');
    fprintf(fid,'end\n');
    fprintf(fid,'contexp(1) = contexp(1) + 1;\n');
    fprintf(fid,'aux = 1;\n');
    fprintf(fid,'while ((aux <= length(vertices)) && (contexp(aux) > size(exponents{aux},1)))\n');
    fprintf(fid,'contexp(aux) = 1;\n');
    fprintf(fid,'contexp(aux+1) = contexp(aux+1) + 1;\n');
    fprintf(fid,'aux = aux + 1;\n');
    fprintf(fid,'end\n');
    fprintf(fid,'end\n');

    fprintf(fid,'return\n\n');
    
    
    
    
    
    if (contZ + contF > 1)
        
        fprintf(fid,'function T = convertaffine(M,bounds)\n\n');
        
        fprintf(fid,'numparam = length(M)-1;\n');
        fprintf(fid,'for cont=1:numparam\n');
        fprintf(fid,'exponents{cont} = [1 0; 0 1];\n');
        fprintf(fid,'end\n');
        
        fprintf(fid,'numelem = 1;\n');
        fprintf(fid,'base = [];\n');
        fprintf(fid,'vertnow = 0;\n');
        fprintf(fid,'for contsimplex=1:numparam\n');
        fprintf(fid,'base = [base (2*ones(1,2)).^(vertnow:vertnow+1)];\n');
        fprintf(fid,'vertnow = vertnow + 2;\n');
        fprintf(fid,'numelem = numelem*size(exponents{contsimplex},1);\n');
        fprintf(fid,'end\n');
        fprintf(fid,'contexp = ones(1,numparam+1);\n');
        fprintf(fid,'for conttotal=1:numelem\n');
        fprintf(fid,'vetexponent = [];\n');
        fprintf(fid,'for contsimplex=1:numparam\n');
        fprintf(fid,'polyexponent{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);\n');
        fprintf(fid,'if (sum(polyexponent{contsimplex}) > 0)\n');
        fprintf(fid,'vetexponent = [vetexponent polyexponent{contsimplex}];\n');
        fprintf(fid,'end\n');
        fprintf(fid,'end\n');
        fprintf(fid,'indhash(conttotal) = sum(base.*vetexponent);\n');
        fprintf(fid,'contexp(1) = contexp(1) + 1;\n');
        fprintf(fid,'aux = 1;\n');
        fprintf(fid,'while ((aux <= numparam) && (contexp(aux) > size(exponents{aux},1)))\n');
        fprintf(fid,'contexp(aux) = 1;\n');
        fprintf(fid,'contexp(aux+1) = contexp(aux+1) + 1;\n');
        fprintf(fid,'aux = aux + 1;\n');
        fprintf(fid,'end\n');
        fprintf(fid,'end\n');
        
        fprintf(fid,'numcoefs = 2^numparam;\n');
        fprintf(fid,'contexp = ones(1,numparam+1);\n');
        fprintf(fid,'T = [];\n');
        fprintf(fid,'for cont=1:numcoefs\n');
        fprintf(fid,'Taux = M{1};\n');
        fprintf(fid,'for contsimplex=1:numparam\n');
        fprintf(fid,'auxexponents{contsimplex} = exponents{contsimplex}(contexp(contsimplex),:);\n');
        fprintf(fid,'Taux = Taux + ((contexp(contsimplex)-1)*bounds(contsimplex,2) + (2-contexp(contsimplex))*bounds(contsimplex,1))*M{contsimplex+1};\n');
        fprintf(fid,'end\n');
        fprintf(fid,'indresul = gethash(auxexponents,indhash,base);\n');
        fprintf(fid,'T{indresul} = Taux;\n');
        
        fprintf(fid,'contexp(1) = contexp(1) + 1;\n');
        fprintf(fid,'aux = 1;\n');
        fprintf(fid,'while ((aux <= numparam) && (contexp(aux) > size(exponents{aux},1)))\n');
        fprintf(fid,'contexp(aux) = 1;\n');
        fprintf(fid,'contexp(aux+1) = contexp(aux+1) + 1;\n');
        fprintf(fid,'aux = aux + 1;\n');
        fprintf(fid,'end\n');
        fprintf(fid,'end\n');
        fprintf(fid,'return\n\n');
        
    end
    
    
    

    fprintf(fid,'\n\n\nfunction index = gethash(exponent,hash,base)\n');
    fprintf(fid,'vetexponent = [];\n');
    fprintf(fid,'for contsimplex=1:(length(exponent)) \n');
    fprintf(fid,'if (sum(exponent{contsimplex}) > 0)\n');
    fprintf(fid,'vetexponent = [vetexponent exponent{contsimplex}];\n');
    fprintf(fid,'end\n');
    fprintf(fid,'end\n');
    fprintf(fid,'index = find(sum(base.*(vetexponent))==hash);\n');
    fprintf(fid,'return\n\n');
    
    
    fprintf(fid,'\n\n\n function exponents = generate_homogeneous_exponents(vertices,degree)\n');
    fprintf(fid,'exponents = [];\n');
    fprintf(fid,'if (vertices == 1)\n');
    fprintf(fid,'exponents = degree;\n');
    fprintf(fid,'return\n');
    fprintf(fid,'end\n');
    fprintf(fid,'if (degree == 0)\n');
    fprintf(fid,'exponents = zeros(1,vertices);\n');
    fprintf(fid,'return\n');
    fprintf(fid,'end\n');
    fprintf(fid,'for cont=degree:-1:0\n');
    fprintf(fid,'[expaux] = generate_homogeneous_exponents(vertices-1,degree-cont);\n');
    fprintf(fid,'exponents = [exponents; cont*ones(size(expaux,1),1) expaux];\n');
    fprintf(fid,'end\n');

    fprintf(fid,'return');
    
  
    fclose(fid);
    fclose(fidtemp);
    
    delete(char(strcat(filename,'.temp')));
    delete('lmifilestemp.mat');
    
    if ~isempty(fidaux)
        fclose(fidaux);
        delete('varaux_temporary_123.temp');
    end
    
end

return







function update_variables(Term,vecpoly)
% Function used to update the informations of all the variables and
% and constants in the generation of the .m file. The informations are
% stored in the file "lmifilestemp", which contains the structure named
% "variables", with fields given by:
% name: Name of the variable
% nummon: Number of monomials
% type: If is variable (V), constant (C), predefined (P), scalar variable
%    (S), scalar constant (K), auxiliar variable (A), affine constant
%    matrix (Z), affine variable matrix (F)
% position: positions [LMI,row,column] where such variable can be found.
% dimension: dimension if it is variable
% varstruct: If it is variable (V), denotes if it is symmetric or full. If
%    it is predefined (P), denotes if it is identity or zero

% The variable "numLMI" informs the number of the current LMI analyzed
% The variable "indexes" informs which indexes contain variables

load lmifilestemp;

% "variables" is a hash table, with key mod(sum(double(variables.name)),1024)+1. In
% case of hit detection, increment the key by one.

% First step: Updating the variable list
for conti=1:size(Term,1)
    LMIrow{conti} = [];
    for contj = conti:size(Term,2)
        LMIcol{contj} = [];
        expr = Term{conti,contj}.opcode{1};
        toprint = [];
        contexpr = 1;
        while (contexpr <= length(expr))
            varnow = [];
            
            while ((contexpr <= length(expr)) && (strcmp(expr(contexpr),'(') || strcmp(expr(contexpr),'-')))
                contexpr = contexpr + 1;
            end
            
            %Name of the variable
            while (~strcmp(expr(contexpr),'#'))
                varnow = strcat(varnow,expr(contexpr));
                contexpr = contexpr + 1;
            end
            vartype = expr(contexpr+1);
            
            key = findhash(varnow,variables);
            
            if (isempty(variables{key}))
                %Insert the new variable
                contpol = 1;
                indpol = 0;
                while (indpol == 0)
                    if (strcmp(vecpoly{contpol}.label,varnow))
                        indpol = contpol;
                    else
                        contpol = contpol + 1;
                    end
                end


                if (issymmetric(vecpoly{indpol}.data(1).value))
                    variables{key}.varstruct = 'symmetric';
                else
                    variables{key}.varstruct = 'full';
                end


                variables{key}.name = varnow;
                if ((vartype == 'F') || (vartype == 'Z'))
                    variables{key}.nummon = vecpoly{indpol}.nummonorig;
                else
                    variables{key}.nummon = length(vecpoly{indpol}.data);
                end
                variables{key}.type = vartype;
                variables{key}.position = [numLMI conti contj];
                if ((vartype == 'V') || (vartype == 'P') || (vartype == 'F') || (vartype == 'D'))
                    variables{key}.dimension{1} = [];
                    variables{key}.dimension{2} = [];
                    indexes = [indexes key];
                end
                if (vartype == 'P')
                    if (sum(sum(vecpoly{indpol}.data(1).value)) == 0)
                        variables{key}.varstruct = 'zero';
                    else
                        variables{key}.varstruct = 'identity';
                    end
                end
            else %It is already inserted, update position
                variables{key}.position = unique([variables{key}.position; numLMI conti contj],'rows');
            end
                
            %Passing the number
            contexpr = contexpr + 2;
            while ((contexpr <= length(expr)) && ~strcmp(expr(contexpr),'+') && ~strcmp(expr(contexpr),'-') && ~strcmp(expr(contexpr),'*') && ~strcmp(expr(contexpr),char(39)) && ~strcmp(expr(contexpr),')') && ~strcmp(expr(contexpr),'('))
                contexpr = contexpr + 1;
            end
            
            %Passing the operations
            while ((contexpr <= length(expr)) && (strcmp(expr(contexpr),'+') || strcmp(expr(contexpr),'-') || strcmp(expr(contexpr),'*') || strcmp(expr(contexpr),char(39)) || strcmp(expr(contexpr),')') || strcmp(expr(contexpr),'(')))
                contexpr = contexpr + 1;
            end
        end
    end
end


%Second step: Updating the dimensions of the variables

changed = false;
contrep = 0; % To test twice if the info on the dimensions were changed
finish = false;
while (~finish)
    contind = 1;
    while (contind <= length(indexes))
        key = indexes(contind);
        contpos = 1;
        while ((contpos <= size(variables{key}.position,1)) && ((isempty(variables{key}.dimension{1})) || (isempty(variables{key}.dimension{2}))))
            %Try to calculate the dimensions of the current variable
            if (variables{key}.position(contpos,1) == numLMI) 
                %If it is in the current LMI
                currow = variables{key}.position(contpos,2);
                curcol = variables{key}.position(contpos,3);
                expr = Term{currow,curcol}.opcode{1};
                
                [row,col,variables] = analyze_expression(expr,variables,LMIrow{currow},LMIcol{curcol});
                %Returns, if possible, the number of rows and columns of
                %this position of the LMI
                
                if ((~isempty(row)) && (isempty(LMIrow{currow})))
                    LMIrow{currow} = row;
                    LMIcol{currow} = row;
                end
                if ((~isempty(col)) && (isempty(LMIcol{curcol})))
                    LMIcol{curcol} = col;
                    LMIrow{curcol} = col;
                end
                
                %If the variable is symmetric or is the identity, equals
                %the number of rows and columns
                if (((variables{key}.type == 'V') && (strcmp(variables{key}.varstruct,'symmetric')))  ||  ((variables{key}.type == 'Z') && (strcmp(variables{key}.varstruct,'symmetric')))  || ((variables{key}.type == 'P') && (strcmp(variables{key}.varstruct,'identity'))))
                    if (~isempty(variables{key}.dimension{1}))
                        variables{key}.dimension{2} = variables{key}.dimension{1};
                    elseif (~isempty(variables{key}.dimension{2}))
                        variables{key}.dimension{1} = variables{key}.dimension{2};
                    end
                end
            end
            contpos = contpos + 1;
        end
        if ((~isempty(variables{key}.dimension{1})) && (~isempty(variables{key}.dimension{2})))
            indexes(contind) = [];
            changed = true;
            contrep = 0;
        else
            contind = contind + 1;
        end
    end
    if (~changed)
        contrep = contrep + 1;
    end
    if ((contrep == 2) || (isempty(indexes)))
        finish = true;
    end
end


numLMI = numLMI + 1;
save lmifilestemp variables numLMI indexes filename
return


function key = findhash(name,variables)
% Function that returns the key of the hash table variables

key = mod(sum(double(name)),1024);
finish = false;
while (~finish)
    if (isempty(variables{key}))
        finish = true;
    elseif (strcmp(name,variables{key}.name))
        finish = true;
    else
        key = key + 1;
        if (key == 1025)
            key = 1;
        end
    end
end

return


function [row,col,variables] = analyze_expression(expr,variables,row,col)
% Function that analyzes an expression and returns, if possible, the
% dimensions of such expression

contexpr = 1;
if (expr(1) == '-')
    contexpr = contexpr + 1;
end

tamstackexpr = 0;
part = [];
%Separate each part between + and - and insert then into stackexpr
while (contexpr <= length(expr))
    if ((expr(contexpr) ~= '+') && (expr(contexpr) ~= '-'))
        part = strcat(part,expr(contexpr));
        if (expr(contexpr) == '(')
            % Obtain the expression on parenthesis
            contexpr = contexpr + 1;
            stackpar = 1;
            while (stackpar > 0)
                if (expr(contexpr) == '(')
                    stackpar = stackpar + 1;
                elseif (expr(contexpr) == ')')
                    stackpar = stackpar - 1;
                end
                part = strcat(part,expr(contexpr));
                contexpr = contexpr + 1;
            end
        end
    else
        tamstackexpr = tamstackexpr + 1;
        stackexpr{tamstackexpr} = part;
        part = [];
    end
    contexpr = contexpr + 1;
end
tamstackexpr = tamstackexpr + 1;
stackexpr{tamstackexpr} = part;
        
%Analyze the number of rows and columns of each factor
transpose = false;
for contstack = 1:length(stackexpr)
    clear rowvar colvar;
    colsafter = [];
    contcolsafter = 1;
    part = stackexpr{contstack};
    contpart = 1;
    numfator = 0;
    prevscalar = false;
    while (contpart < length(part))
        numfator = numfator + 1;
        rowvar{numfator} = [];  colvar{numfator} = [];
        fator = [];
        if (part(contpart) == '(')
            % Obtain the expression on parenthesis and analyze recursively
            contpart = contpart + 1;
            stackpar = 1;
            while (stackpar > 0)
                if (part(contpart) == '(')
                    stackpar = stackpar + 1;
                elseif (part(contpart) == ')')
                    stackpar = stackpar - 1;
                end

                if (stackpar > 0)
                    fator = strcat(fator,part(contpart));
                end
                contpart = contpart + 1;
            end
            [rowvar{numfator},colvar{numfator},variables] = analyze_expression(fator,variables,row,col);
        else
            %Name of the variable
            while (~strcmp(part(contpart),'#'))
                fator = strcat(fator,part(contpart));
                contpart = contpart + 1;
            end
            vartype = part(contpart+1);

            %Passing the number
            contpart = contpart + 2;
            while ((contpart <= length(part)) && ~strcmp(part(contpart),'*') && ~strcmp(part(contpart),char(39)) && ~strcmp(part(contpart),')') && ~strcmp(part(contpart),'('))
                contpart = contpart + 1;
            end

            %Passing the operations
            transpose = false;
            while ((contpart <= length(part)) &&  (strcmp(part(contpart),'*') || strcmp(part(contpart),char(39)) || strcmp(part(contpart),')') || strcmp(part(contpart),'(')))
                if (part(contpart) == char(39))
                    transpose = true;
                end
                contpart = contpart + 1;
            end

            if ((vartype == 'S') || (vartype == 'K'))
                if (numfator == 1)
                    rowvar{numfator} = row;
                    colvar{numfator} = row;
                else
                    rowvar{numfator} = colvar{numfator-1};
                    colvar{numfator} = colvar{numfator-1};
                end
                if (isempty(rowvar{numfator}))
                    prevscalar = true;
                end
            elseif ((vartype == 'C') || (vartype == 'Z'))
                if (~transpose)
                    rowvar{numfator} = strcat('row',fator);
                    colvar{numfator} = strcat('col',fator);
                else
                    rowvar{numfator} = strcat('col',fator);
                    colvar{numfator} = strcat('row',fator);
                end
                if (numfator == 1)
                    row = rowvar{numfator};
                end
                if (prevscalar)
                    rowvar{numfator-1} = rowvar{numfator};
                    colvar{numfator-1} = colvar{numfator};
                end
                prevscalar = false;
            else
                key = findhash(fator,variables);
                if (~transpose)
                    rowvar{numfator} = variables{key}.dimension{1};
                    colvar{numfator} = variables{key}.dimension{2};
                else
                    rowvar{numfator} = variables{key}.dimension{2};
                    colvar{numfator} = variables{key}.dimension{1};
                end
                if ((numfator == 1) && (~isempty(rowvar{numfator})))
                    row = rowvar{numfator};
                end
                if (prevscalar && (~isempty(rowvar{numfator})))
                    rowvar{numfator-1} = rowvar{numfator};
                    colvar{numfator-1} = colvar{numfator};
                end
                prevscalar = false;
                
                %Analyze the row
                changed = false;
                if (isempty(rowvar{numfator}))
                    if ((numfator > 1) && (~isempty(colvar{numfator-1})))
                        rowvar{numfator} = colvar{numfator - 1};
                        changed = true;
                    elseif ((numfator == 1) && (~isempty(row)))
                        rowvar{numfator} = row;
                        changed = true;
                    end
                end
                if (changed)
                    if (~transpose)
                        variables{key}.dimension{1} = rowvar{numfator};
                    else
                        variables{key}.dimension{2} = rowvar{numfator};
                    end
                end
                
                %Store the columns to be analyzed
                if (isempty(colvar{numfator}))
                    colsafter = [colsafter numfator];
                    fatorsafter{contcolsafter} = fator;
                    transposeafter(contcolsafter) = transpose;
                    contcolsafter = contcolsafter + 1;
                end
            end
        end
    end
    %Analyze the columns of the variables
    totalfators = numfator;
    if (prevscalar)
        %The last variable is a scalar
        rowvar{totalfators} = col;
        colvar{totalfators} = col;
    end
    for contcol = 1:length(colsafter)
        numfator = colsafter(contcol);
        fator = fatorsafter{contcol};
        transpose = transposeafter(contcol);

        changed = false;
        if ((numfator < totalfators) && (~isempty(colvar{numfator+1})))
            colvar{numfator} = rowvar{numfator + 1};
            changed = true;
        elseif ((numfator == totalfators) && (~isempty(col)))
            colvar{numfator} = col;
            changed = true;
        end
        if (changed)
            key = findhash(fator,variables);
            if (~transpose)
                variables{key}.dimension{2} = colvar{numfator};
            else
                variables{key}.dimension{1} = colvar{numfator};
            end
        end
    end
    %Passed by all the factors in this part of the stack
    if (~isempty(colvar{end}))
        if (~transpose)
            col = colvar{end};
        else
            col = rowvar{end};
        end
    end
end

return