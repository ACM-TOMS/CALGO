function [m,D,H]=multiplicity(f,variables,zero,options)
%  MULTIPLICITY --> Computing the multiplicity and multiplicity structure
%  of a system of nonlinear equations at an isolated zero.
%
%  <Synopsis>
%          m = multiplicity(f, variables, zero)
%          m = multiplicity(f, variables, zero, options)
%  [m, D, H] = multiplicity(f, variables, zero)
%  [m, D, H] = multiplicity(f, variables, zero, options)
%
%  <Input Parameters>
%    1. f         --> a cell array containing the system of nonlinear
%                     equations as strings, e.g.,
%                         >>   f = { 'x^2 + sin(y) -1',  'x-cos(y)' };
%
%    2. variables --> a cell array containing the unknown variables as
%                     strings, e.g.,
%                         >>   variables = { 'x',  'y' };
%
%    3. zero      --> a vector containing numerical zero (root) of f, e.g.,
%                         >>   zero = [1, 0];
%
%    4. options   --> an optional parameter which includes the configuration
%                     settings:
%                     Display: Set to 2 to have all output (the dual space
%                              and Hilbert function) printed to the screen,
%                              and set to 1 to have depth and Hilbert
%                              function printed to the screen. Otherwise
%                              set the default value 0;
%                         Tol: The threshold for numerical rank-revealing.
%                              Singular values above Tol are counted as
%                              nonzero. The default value is 1e-8;
%                     EqnType: The equation type for MULTIPLICITY,
%                              polynomial system ('Poly') or nonlinear
%                              functions ('Nonlinear'). The default value is
%                              'Nonlinear'. By setting the value to 'Poly',
%                              MULTIPLICITY will transfer the polynomial
%                              system to the matrix representation
%                              internally and speed up the computation.
%                      MaxInt: Maximum multiplicity allowed in the recursive
%                              computation. If a zero is not isolated, the
%                              multiplicity will be infinity.  The code can
%                              be used for identifying a nonisolated zero by
%                              setting MaxInt to a known upper bound (e.g.
%                              the Bezout number). The default value is 1000.
%
%                     All the configuration settings may be changed by using
%                     optset function, i.e.,
%                         >>   options = optset('para', value);
%                     para could be any configuration setting, value is set
%                     to para. See OPTSET for details. Any configuration
%                     setting that is not changed will be set to its default
%                     value.
%
%  <Output Parameters>
%    1. m        --> the multiplicity of f at the zero;
%
%    2. D        --> a basis for the dual space as a cell array with
%                    each component being a matrix in the Matlab form of
%                                   D{i} = [c_1, j_1;
%                                           c_2, j_2;
%                                             ...
%                                           c_n, j_n ];
%                    representing a differential functional
%                            D_i = c_1*d_{j_1} + ... + c_n*d_{j_n}
%                    where d_{j_i}'s are differential monomial functionals
%                    (e.g. For a system of equations with variables {x,y,z}
%                     at the zero x=a, y=b, z=c, the functional d_{[i,j,k]}
%                     applied to any function g is the value of the partial
%                     derivative
%
%                                                         i+j+k
%                                             1          d
%                         d_{[i,j,k]}(g) = -------- * ----------- g(a,b,c)
%                                          i! j! k!     i   j   k
%                                                     dx  dy  dz
%
%                     The dual space is the vector space consists of such
%                     differential functionals that vanish on the nonlinear
%                     system while satisfying the so-called closedness
%                     condition);
%
%    3. H        --> values of the Hilbert function in a vector.
%
%  <Examples>
%   Consider the nonlinear system
%
%               sin(x)*cos(y)-x       = 0
%               sin(y)*sin(x)^2 - y^2 = 0
%
%   at the zero (x, y) = (0, 0), the multiplicity can be computed by the
%   following statements:
%
%     >> f = {'sin(x)*cos(y) - x', 'sin(y)*sin(x)^2 - y^2'};
%     >> variables = {'x', 'y'};
%     >> zero = [0, 0];
%     >> m = multiplicity(f, variables, zero)
%
%   To create an options structure with Tol = 1e-10, MaxInt = 100:
%     >> options = optset('Tol', 1e-10, 'MaxInt', 100);
%     >>  m = multiplicity(f, variables, zero, options)
%
%  <Algorithm>
%   This code implements a modified closedness subspace method for
%   multiplicity identification with a newly developed equation-by-equation
%   strategy for improving efficiency.
%
%  <References>
%  [1] An algorithm and software for computing multiplicity structures at
%      zeros of nonlinear systems, W. Hao, A. J. Sommesse and Z. Zeng
%  [2] The closedness subspace method for computing the multiplicity
%      structure of a polynomial system, Z. Zeng.
%
%  See also optset, cell

%Set the default configuration
if nargin==3
    options=defaultstructure;
end

if strcmp(options.EqnType,'Poly')
    %the input f is polynomial, setup f and variables
    matrixf=nonlinear2matrix(f,variables);
    [m,D,H]=multiplicity_poly(matrixf,zero,options.Tol,...
        options.MaxInt,options.Display);
else
    %the input f is nonlinear function, setup f, variables and x
    [m,D,H]=multiplicity_nonlinear(f,variables,zero,options.Tol,...
        options.MaxInt,options.Display);
end
end

function [mul,D,H]=multiplicity_poly(f,x,tol,MaxInt,dispflag)
% multiplicity_poly(f,x): determine the multiplicity structrue of
% a polynomial system f at an isolated  zero.
%
% Input:
%   f =  an arrary, which represent the polynomial system
%       {f_1,f_2,...,f_t}
%   x =  an isolated zero for system f
%   tol = the threshold tolerance of zero,default value is 1e-8
%   dispflag = 0  no print (default)
%              1  print the depth and multiplicity
%              2  print the basis for each dual space
% Output:
%  mul = the multiplicity
%   D  = the dual space of f at x
%   H  = the hilbert function
%

if length(x)~=size(f,2)-1
    error('The number of variables in x and f are not same!')
end
% Initalize the basis for dual subspace D and closedness subspace Ca
Sa0=1;
Ma1=[zeros(size(x));eye(size(x,2))];% Initalize closedness support
m=max(max(Ma1))+1;
IndexMa=Ma1*m.^(0:size(Ma1,2)-1)'; % Indices of monomial in support
numMa1=1;                  % Number of previous monomial in support
Ca0=[1 zeros(size(x))];    % Previous closedness subspace
N0=1;                      % Compute null space for Wa
alpha=0;                   % Depth for dual spaces
DI{alpha+1}=1;             % Initalize dual subspace
df_count=1;
if dispflag
    disp('   Depth     Hilbert function')
end
mul=0;
while(1)
    %
    % The main loop of the process
    %
    Ca1=generateCa(Ma1,DI,Ca0,Sa0,m,IndexMa,tol);
    % Get updated closedness subspace Ca by the algorithm in Sec. 6
    % Print to screen if Display option is 1 or 2
    if dispflag
        disp(['      ',num2str(alpha),'                 ',...
            num2str(size(DI{alpha+1},1))])
    end
    mul=size(DI{alpha+1},1)+mul;  % Update multiplicity count
    if mul>MaxInt  % exit if multiplicity is larger than MaxInt
        break
    end
    alpha=alpha+1;% Depth increases by 1
    if size(Ca1,1)==size(Ca0,1)
        % the process is completed when dim(Ca0)=dim(Ca1), break the loop
        break
    end
    %
    for i=numMa1:size(Ma1,1)
        %
        % loop for updating df
        %
        if df_count==1
            df = ployeval(dfpoly(f,Ma1(i,:)),x);
        else
            df(:,df_count)=ployeval(dfpoly(f,Ma1(i,:)),x);
        end
        df_count=df_count+1;
    end
    %
    % Calculate the null space of Wa
    %
    Wa=df*Ca1';   % Generate Wa
    N0(size(N0,1)+1:size(Ca1,1)-size(Ca0,1)+size(N0,1),:)=...
        zeros(size(Ca1,1)-size(Ca0,1),size(N0,2));
    Waorth=mynull([ N0'; Wa],tol);  % Calculate null space by mynull
    N=[N0 Waorth];                  % Update null space
    N(abs(N)<tol)=0;                % Clear tiny entries in N
    if size(N,2)==size(N0,2)
        % Null space stops expanding, the process is completed
        break
    else
        % otherwise, update N0 for next step
        N0=N;
    end
    Da=N'*Ca1;           % Generate new Da
    Da(abs(Da)<tol)=0;   % Clear tiny entries in Da
    index=Da~=0;
    index=find(sum(index)~=0);
    Sa1=Sa0;
    Sa0=union(Sa1,index);
    [Sa1,index1]=intersect(Sa0,Sa1);
    % Update the dual subspace DI
    tmp=1;
    if length(Sa1)<length(Sa0)
        for i=1:alpha
            DItmp=zeros(size(DI{i},1),length(Sa0));
            DItmp(:,index1)=DI{i};
            DI{i}=DItmp;
            tmp=tmp+size(DI{i},1);
        end
        DI{alpha+1}=Da(tmp:end,index);
    else
        break
    end
    % Update monomial support
    numMa1=size(Ma1,1)+1;
    % Prepare for the next iteration, update Ma1 and Ca0
    Ma1=expandMa(Ma1,Sa0,Sa1,size(x,2),m,IndexMa);
    m=max(max(Ma1))+1;
    IndexMa=Ma1*m.^(0:size(Ma1,2)-1)';
    % Save the current closedness subspace
    Ca0=[Ca1 zeros(size(Ca1,1),size(Ma1,1)-size(Ca1,2))];
end
% process output
if dispflag
    disp('---------------------------------------------')
    disp(['Depth:',num2str(alpha-1),'    Multiplicity:',num2str(mul)]);
end
%output the detail information if needed
if mul>MaxInt
    disp(['The multiplicity exceeds MaxInt.',...
        'The zero might be a nonisolated solution!',...
        ' Or set ''MaxInt'' be a largger number!'])
end
Sa0=Ma1(Sa0,:);
if dispflag>1
    disp('The dual space is:')
    muloutput(DI,Sa0)
end
%output DI HF
count=1;
H=zeros(1,length(DI));
D=cell(1,length(DI)*size(DI{end},1));
for i=1:length(DI)
    for j=1:size(DI{i},1)
        tmp=find(DI{i}(j,:)~=0);
        D{count}=[DI{i}(j,tmp)' Sa0(tmp,:)];
        count=count+1;
    end
    H(i)=size(DI{i},1);
end
D=D(1:count-1);
end

function [mul,D,H]=multiplicity_nonlinear(f,var,x,tol,MaxInt,dispflag)
% multiplicity_nonlinear(f,x): determine the multiplicity structrue of a
% nonlinear system f at an isolated  zero.
%
% Input:
%   f   = a cell structure, which represent the nonlinear system
%   var = the variable name list
%   x   = an isolated zero for system f.
%   tol = the threshold tolerance of zero,default value is 1e-8
%   dispflag = 0  no print (default)
%              1  print the depth and multiplicity
%              2  print the basis for each dual space
% Output:
%   mul = the multiplicity
%   D   = the dual space of f at x
%   H   = the hilbert function

if length(x)~=length(f)
    error('The number of variables in x and f are not same!')
end
%Evaluate the variables
for i=1:length(x)
    if i==1
        if iscell(var)
            ff=[var{i},'=', num2str(x(i),16),';'];
        else
            ff=[char(var(i)),'=', num2str(x(i),16),';'];
        end
        f_count=length(ff);
    else
        if iscell(var)
            ftmp=[var{i},'=', num2str(x(i),16),';'];
        else
            ftmp=[char(var(i)),'=', num2str(x(i),16),';'];
        end
        ff(f_count+1:length(ftmp)+f_count) = ftmp ;
        f_count=length(ftmp)+f_count;
    end
end
% Initalize the basis for dual subspace D and closedness subspace Ca
Sa0=1;
Ma1=[zeros(size(x));eye(size(x,2))];% Initalize closedness support
m=max(max(Ma1))+1;
IndexMa=Ma1*m.^(0:size(Ma1,2)-1)'; % Indices of monomial in support
numMa1=1;                  % Number of previous monomial in support
Ca0=[1 zeros(size(x))];    % Previous closedness subspace
N0=1;                      % Compute null space for Wa
alpha=0;                   % Depth for dual spaces
DI{alpha+1}=1;             % Initalize dual subspace
df_count=1;
if dispflag
    disp('   Depth     Hilbert function')
end
mul=0;
while(1)
    %
    % The main loop of the process
    %
    Ca1=generateCa(Ma1,DI,Ca0,Sa0,m,IndexMa,tol);
    % Get updated closedness subspace Ca by the algorithm in Sec. 6
    % Print to screen if Display option is 1 or 2
    if dispflag
        disp(['      ',num2str(alpha),'                 ',...
            num2str(size(DI{alpha+1},1))])
    end
    mul=size(DI{alpha+1},1)+mul;  % Update multiplicity count
    if mul>MaxInt    % exit if multiplicity is larger than MaxInt
        break
    end
    alpha=alpha+1;%Depth increases by 1
    if size(Ca1,1)==size(Ca0,1)
        % The process is completed when dim(Ca0)=dim(Ca1), break the loop
        break
    end
    for i=numMa1:size(Ma1,1)
        %
        % loop for updating df
        %
        if df_count==1
            df = diff_f(f,var,Ma1(i,:),ff);
        else
            df(:,df_count)=diff_f(f,var,Ma1(i,:),ff);
        end
        df_count=df_count+1;
    end
    %
    % Calculate the null space of Wa
    %
    Wa=df*Ca1';  % Generate Wa
    N0(size(N0,1)+1:size(Ca1,1)-size(Ca0,1)+size(N0,1),:)=...
        zeros(size(Ca1,1)-size(Ca0,1),size(N0,2));
    Waorth=mynull([ N0'; Wa],tol);  % Calculate null space by mynull
    N=[N0 Waorth];                  % Update null space
    N(abs(N)<tol)=0;                % Clear tiny entries in N
    if size(N,2)==size(N0,2)
        % Null space stops expanding, the process is completed
        break
    else
        % Otherwise, update N0 for next step
        N0=N;
    end
    Da=N'*Ca1;           % Generate new Da
    Da(abs(Da)<tol)=0;   % Clear tiny entries in Da
    index=Da~=0;
    index=find(sum(index)~=0);
    Sa1=Sa0;
    Sa0=union(Sa1,index);
    [Sa1,index1]=intersect(Sa0,Sa1);
    % Update the dual subspace DI
    tmp=1;
    if length(Sa1)<length(Sa0)
        for i=1:alpha
            DItmp=zeros(size(DI{i},1),length(Sa0));
            DItmp(:,index1)=DI{i};
            DI{i}=DItmp;
            tmp=tmp+size(DI{i},1);
        end
        DI{alpha+1}=Da(tmp:end,index);
    else
        break
    end
    % Update monomial support
    numMa1=size(Ma1,1)+1;
    % Prepare for the next iteration, update Ma1 and Ca0
    Ma1=expandMa(Ma1,Sa0,Sa1,size(x,2),m,IndexMa);
    m=max(max(Ma1))+1;
    IndexMa=Ma1*m.^(0:size(Ma1,2)-1)';
    Ca0=[Ca1 zeros(size(Ca1,1),size(Ma1,1)-size(Ca1,2))];
    % Save the current closedness subspace
end
% Process output
if dispflag
    disp('---------------------------------------------')
    disp(['Depth:',num2str(alpha-1),'    Multiplicity:',num2str(mul)]);
end
%output the detail information if needed
if mul>MaxInt
    disp(['The multiplicity exceeds MaxInt.',...
        'The zero might be a nonisolated solution!',...
        'Or set ''MaxInt'' be a largger number!'])
end
Sa0=Ma1(Sa0,:);
if dispflag>1
    disp('The dual space is:')
    muloutput(DI,Sa0)
end
%output DI HF
count=1;
H=zeros(1,length(DI));
D=cell(1,length(DI)*size(DI{i},1));
for i=1:length(DI)
    for j=1:size(DI{i},1)
        tmp=find(DI{i}(j,:)~=0);
        D{count}=[DI{i}(j,tmp)' Sa0(tmp,:)];
        count=count+1;
    end
    H(i)=size(DI{i},1);
end
D=D(1:count-1);
end

function matrixf=nonlinear2matrix(f,var)
%
% Convert a polynomial system to coefficient matrix
%
eof=zeros(1,length(var)+1);
eof(2)=-1;
% Set end of equation.
matrixf_count=1;
% Matrix format
if ~iscell(f)
    f1=cell(1,length(f));
    for ff=1:length(f)
        f1{ff}=char(f(ff));
    end
    f=f1;
    clear f1;
end
if ~iscell(var)
    var1=cell(1,length(var));
    for ff=1:length(var)
        var1{ff}=char(var(ff));
    end
    var=var1;
    clear var1;
end
for ff=1:length(f)
    f{ff}=char(simplify(sym(f{ff})));
    findex= f{ff}==' ';
    f{ff}(findex)=[];
    findex=find(f{ff}=='+');
    findex=union(findex,find(f{ff}=='-'));
    findex=sort(union(union(findex,1),length(f{ff})+1));
    %find the index for each monomial
    for i=1:length(findex)-1
        %pick one monomial
        monial=f{ff}(findex(i):findex(i+1)-1);
        monialindex=[0 find(monial=='*') length(monial)+1];
        %inital the coefficient
        if monial(1)=='-'
            coef=str2double([monial(1) '1']);
        else
            coef=1;
        end
        %one row of matrix
        tmpmf=zeros(1,length(var)+1);
        for j=1:length(monialindex)-1
            %get each term of the monomial
            monterm=monial(monialindex(j)+1:monialindex(j+1)-1);
            iscoef=1;
            for k=1:length(var)
                varindex=strfind(monterm,var{k});
                %find the position
                if length(varindex)==1
                    iscoef=0;
                    break
                end
            end
            if iscoef
                %evulate the coefficient
                coef=str2double(monterm);
            else
                %get the power
                if length(monterm)>=varindex+length(var{k})
                    if monterm(varindex+length(var{k}))=='^'
                        tmpmf(k+1) = str2double(monterm(varindex+...
                            length(var{k})+1:end));
                    else
                        tmpmf(k+1)=1;
                    end
                else
                    tmpmf(k+1)=1;
                end
            end
        end
        tmpmf(1)=coef;
        if matrixf_count==1
            matrixf=tmpmf;
            matrixf_count=matrixf_count+size(tmpmf,1);
        else
            matrixf(matrixf_count:matrixf_count+size(tmpmf,1)-1,:)=tmpmf;
            matrixf_count=matrixf_count+size(tmpmf,1);
        end
    end
    matrixf(matrixf_count,:)=eof;
    matrixf_count=matrixf_count+1;
end
end

function df=dfpoly(f,dx)
% df=dfpoly(f,dx)  computes derivative for f, which is a polynomial system.
% dx is a array representing  \partial_{j_1,...,j_s}
% df is the polynomial system after derviatives
%
eof=zeros(1,length(dx)+1);
eof(2)=-1;%ending polynomial equation sign.
df=eof;
tmp=1;%The number of polynomial equations
for i=1:size(f,1)
    %check whether derivative of fi is 0, if it is 0, dfi=zero, else
    %df=[df,eof],which is notation for end of poly equation
    if f(i,:)==eof
        if tmp==1
            df=zeros(size(eof));
            tmp=tmp+1;
        else
            if df(tmp-1,:)==eof
                df(tmp,:)=zeros(size(eof));
                tmp=tmp+1;
            end
        end
        df(tmp,:)=eof;
        tmp=tmp+1;
        continue
    end
    %loop: derivatives for every xi
    df(tmp,2:end)=f(i,2:end)-dx;
    if min(df(tmp,2:end))<0
        continue
    else
        df(tmp,1)=prod([f(i,1) CombNM(f(i,2:end)-dx,f(i,2:end))]);
        tmp=tmp+1;
    end
end
df(:,1)=df(:,1)/prod(factorial(dx));
%factorial mulitplication for every j_i
end

function df=diff_f(f,var,Ma,ff)
%Get the derivative of function f to Ma.
df=sym(zeros(length(f),1));
for i=1:length(f)
    if iscell(f)
        df(i)=diff(sym(f{i}),sym(var{1}),Ma(1));
    else
        if iscell(var)
            df(i)=diff(f(i),sym(var{1}),Ma(1));
        else
            df(i)=diff(f(i),var(1),Ma(1));
        end
    end
end
for i=2:length(var)
    if iscell(var)
        df=diff(df,sym(var{i}),Ma(i));
    else
        df=diff(df,var(i),Ma(i));
    end
end
%Get the value of x
eval(ff)
%Evaluate df at x
df=eval(df)/prod(factorial(Ma));
end

function Ma=expandMa(Ma1,Sa1,Sa2,s,m,IndexMa)
% Ma=expandMa(Ma1,Sa1,Sa2,s) set up the closedness support M^{\alpha} of
% order \alpha from the dual support S^{\alpha-1} of order \alpha-1
%
% Input:
%   Ma1 : Closedness support M^{\alpha-1} of order \alpha-1
%   Sa1 : Dual support S^{\alpha-1} of order \alpha-1
%   Sa2 : Dual support S^{\alpha-2} of order \alpha-2
%     s : Dimensions of variables
%
%Output:
%    Ma : Closedness support M^{\alpha} of order \alpha
%

Sa12=Sa1;%Generate the S^{\alpha-1}\S^{\alpha-2}
[~,Id]=intersect(Sa1,Sa2);
Sa12(Id)=[];
if isempty(Sa12)
    Ma=Ma1;%If Sa12=empty, M^{\alpha}=M^{\alpha-1} return
    return
end
T=[];%Initaialize the temporary set
for i=1:s
    %Step 1 of alogrithm (21) create the set of possible monomials
    T=unique([T;psaidelta(Ma1(Sa12,:),i)],'rows');
end
count=1:size(T,1);
for i=1:size(T,1)
    %exclude monomials not satisfying the closedness condition
    flag=1; % Whether it satisify closedness condition
    for j=1:s
        tmp=feidelta(T(i,:),j,m);
        if tmp~=-1
            tmp=find(IndexMa(Sa1)==tmp, 1);
            if isempty(tmp)
                flag=0;
                break
            end
        end
    end
    if flag==0
        %not satisfying the closedness condition, remove it from T
        count(i)=0;
    end
end
Ma=[Ma1;T(count~=0,:)];%return Ma
end

function A=feidelta(A,i,m)
%Linear anti-differentiation operator based on the closedness condition
if A(i)==0
    A=-1;
    return
else
    A(i)=A(i)-1;
    A=A*m.^(0:size(A,2)-1)';
    return
end
end

function Ca=generateCa(Ma,Da1,Ca1,Sa0,mn,IndexMa,tol)
% Ca = generateCa(Ma,Da1,Sa1,Ca1)
% Generating new closedness subspace C^{\alpha} according to section 4's
% algorithm
%
%Input:
%    Ma is closedness support of order \alpha
%   Da1 is the dual space of order \alpha-1
%   Sa1 is the dual support of order \alpha-1
%   Ca1 is the closedness subspace of order \alpha-1
%
%Output:
%    Ca is the closedness subspace of order \alpha
%

Ma(1,:)=[];%Removing \partial_{0 ,....,0}
Ca=Ca1;
Ca(:,1)=[];
%Deleting coefficient corresponding to b
Ca(1,:)=[];
%Compute the number of basises for D^{\alpha}
for ii=1:size(Da1,2)
    if ii==1
        D=Da1{ii};
    else
        D(size(D,1)+1:size(D,1)+size(Da1{ii},1),:)=Da1{ii};
    end
end
D=D';
[m,s]=size(Ma);
n=size(D,2);
for delta=1:s
    tmpr_count=1;
    A=zeros(m,n);
    Aindex=1;
    B=zeros(size(D,1),n);
    %Compute the matrix A,B which is mentioned in Section 4
    Bindex=1;
    indexZ_count=1;
    for j=1:m
        tmp=feidelta(Ma(j,:),delta,mn);
        if tmp~=-1%Check feidelta(Ma(j,:),delta)==0?
            tmp=find(IndexMa(Sa0)==tmp);
            %Find feidelta(Ma(j,:),delta)'s position in S^{\alpha-1}
            A(Aindex,:)=D(tmp,:);%no, evaluation directly
            Aindex=Aindex+1;
            if tmpr_count==1
                tmpr=tmp;
                tmpr_count=tmpr_count+length(tmp);
            else
                tmpr(tmpr_count:tmpr_count+length(tmp)-1)=tmp;
                tmpr_count=tmpr_count+length(tmp);
            end
            %Record the basis in S^{\alpha-1} which are checked
            if indexZ_count==1
                indexZ=j;
                indexZ_count=indexZ_count+1;
            else
                indexZ(indexZ_count)=j;
                indexZ_count=indexZ_count+1;
            end
        end
    end
    if size(tmpr,2)~=size(D,1)%Is all basis in S^{\alpha-1} checked?
        for ii=1:size(D,1)
            % if no, put the coefficients of the basis which is not checked
            % into B
            if isempty(find(tmpr==ii, 1))
                B(Bindex,:)=D(ii,:);
                Bindex=Bindex+1;
            end
        end
    end
    A=A(1:Aindex-1,:);
    B=B(1:Bindex-1,:);
    %Using equation-by-equation method
    if delta==1
        %generate Mim1 directly for delta=1
        Mim1=A*mynull(B,tol);
        indexZi=indexZ;
    else
        %Find the intersection
        [~,indexD1,indexD2]=intersect(indexZi,indexZ);
        D1=Mim1(indexD1,:);
        Mim1(indexD1,:)=[];
        M=A*mynull(B,tol);
        D2=M(indexD2,:);
        %Compute Mi
        Mi=mynull([D1 -D2],tol);
        %Update Mim1
        Mim11=zeros(size(Mim1,1)+size(M,1),size(Mi,2));
        Mim11(1:size(Mim1,1),:)=Mim1*Mi(1:size(Mim1,2),:);
        Mim11(size(Mim1,1)+1:end,:)=M*Mi(1+size(Mim1,2):end,:);
        Mim1=Mim11;
        %Union the index for Zi
        indexZi(indexD1)=[];
        indexZi(length(indexZi)+1:length(indexZi)+length(indexZ))=...
            indexZ;
    end
end
%Rearrange Mim1 follow the order of Zi
[~,I]=sort(indexZi);
Mim1=Mim1(I,:);
Ca1orth=Mim1*mynull(Ca*Mim1,tol);
%The orthogonal space of C^{\alpha-1}
Ca=[Ca1;[zeros(1,size(Ca1orth,2));Ca1orth]'];
%Inserting C^{\alpha-1} into C^{\alpha}
Ca(abs(Ca)<tol)=0;
[~,Caind]=intersect(Ca,zeros(1,size(Ca,2)),'rows');
Ca(Caind,:)=[];
end

function muloutput(DI,DIbasis)
basis=sym(zeros(size(DIbasis,1),1));
for i=1:size(DIbasis,1)
    tmp=num2str(DIbasis(i,:),'%d');
    tmp(tmp==' ')=[];%remove the space
    basis(i)=sym(['d' tmp]);
end
for i=1:length(DI)
    Da=DI{i}*basis;
    tmph=['   D',num2str(i-1),':'];
    for j=1:size(Da,1)
        if j==1
            tmp=[tmph,char(vpa(Da(j,:),16))];
            tmph=ones(1,length(tmph))*' ';
        else
            tmp=[tmph,char(vpa(Da(j,:),16))];
        end
        disp(tmp)
    end
end
end

function N=mynull(A,tol)
%
% calculate numerical null space of a matrix within tol
%
if nargin==1
    tol=1e-8;
end
[m,n]=size(A);
if m>n
    R=qr(A,0);
    [~,S,V] = svd(triu(R(1:n,1:n)),0);
    r = sum(diag(S) > tol);
    N = V(:,r+1:n);
else
    [~,S,V] = svd(A,0);
    r = sum(diag(S) > tol);
    N = V(:,r+1:n);
end

end

function y=ployeval(f,x)
%Evaluate of polynomial system at x
tmp=0;
y_count=1;
for i=1:size(f,1)
    if f(i,2)==-1
        if y_count==1
            y=tmp;
            y_count=y_count+1;
        else
            y(y_count)=tmp;
            y_count=y_count+1;
        end
        tmp=0;
        continue
    end
    tmp=tmp+prod([f(i,1) x].^[1 f(i,2:end)]);
end
y=y';
end

function A1=psaidelta(A,i)
%psaidelta(A,i) linear differentiation operator
%                      \partial_{j_1',...,j_s'}
%Where j_t'=j_t for j \in{1,..,s}\{i}, j_i'=j_i
%Input:
% Matrix A:  each rows represent [j_1',...,j_s']
% i: for i-th differentiation operator.
%
A1=A+[zeros(size(A,1),i-1),...
    ones(size(A,1),1),zeros(size(A,1),size(A,2)-i)];
end

function C=CombNM(N,M)
C=zeros(1,length(N));
for i=1:length(N)
    if M(i)==0 && N(i)==0
        C(i)=1;
    else
        C(i)=prod(N(i)+1:M(i));
    end
end
end

function options=defaultstructure()
options.EqnType='Nonlinear';
options.Display=0;
options.Tol=1e-8;
options.MaxInt=1000;
end