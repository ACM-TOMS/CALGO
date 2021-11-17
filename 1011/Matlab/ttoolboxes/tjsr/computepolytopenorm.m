function [ normval, iteration ] = computepolytopenorm( pts, VV, algorithm, numcore, epsequal, verbose, solver );
% [ normval, iteration ] = computepolytopenorm ( pts, VV, algorithm, [numcore], [epsequal], [verbose], [solver] );
% Computes the Minkowski-norm of a set of points with respect to the polytope co_* VV, where co_* is some convex hull defined by <algorithm>
%
% The function contains one method using MATLAB linprog, and one method using GUROBI solver
% To change from one to the other, one has to comment out/in some lines in the first function
%
% Input:
%   pts             the pts to be tested
%   VV              the vertices of the polytope. Each vertex is one column of the matrix
%   algorithm       which algorithm/convex hull to use
%                       0: cone (co_-)
%                       1: balanced polytope (minkowski) (co_s)
%                       2: abs-polytope (complex) (absco) (not implemented yet)
%   epsequal        (optional) only has an effect on the text-output, determines, when a vertex is considered to be outside
%   verbose         (optional) verbose level
%   numcore         (optional) number of threads used for computation. If numcore==1, then computation is done in main thread (for gurobi solver)
%   
% Output:
%   normval         The norms of the points as a vector
%   iteration       The number of iterations the solver has done
%
% Info:
%   a) This function is optimized to compute the norms of many points wrt to the same polytope
%      If the norm of only one point shall be computed, this function is very slow.
%
% 	b) The type of polytope which is defined by the vertices in <VV>, is defined by <algorithm>
%   c) epsequal only has an effect on the text-output of this function. For the computation a fixed accuracy is used, depending on the algorithm.
%   d) The functions text-output does not end with a newline-character.
%
% E.g.: [normval,iteration]=computepolytopenorm([1 1; 1 0], [1 0; 0 1], 1, 0, 2)
%   
% Written by tommsch, 2018

% Make some basic tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( nargin<=6 || isempty(solver) ); 
    solver = 'a'; end;
if( nargin<=5 || isempty(verbose) ); 
    verbose = 1; end;
if( nargin<=4 || isempty(epsequal) ); 
    epsequal = 0; end;
if( nargin<=3 || isempty(numcore) || isequal(numcore,0) ); 
    numcore = feature( 'numcores' ); end;

if( isempty(pts) ); 
    normval = zeros( 1, 0 ); 
    iteration = 0; 
    return; end;
if(isempty(VV)); 
    normval = inf( 1, size(pts,2) ); 
    iteration = 0; 
    return; end;
if( algorithm~=TJSR_CONEFUNCT && rank(VV)<size(VV,1) );
    normval = inf( 1, size(pts,2) ); 
    iteration = 0; 
    return; end;

% Start the computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch solver
    case {'gurobi','g'};
        [normval,iteration] = computepolytopenorm_gurobi( pts, VV, algorithm, numcore, epsequal, verbose ); %Use Gurobi
    case {'matlab','m'};
        [normval,iteration] = computepolytopenorm_matlab( pts, VV, algorithm, numcore, epsequal, verbose );
    case {'auto','a'}; %automatic
        if( isequal(exist('gurobi_setup','file'),2) && isequal(exist('gurobi','file'),2) );
            try;
                [normval,iteration] = computepolytopenorm_gurobi( pts, VV, algorithm, numcore, epsequal, verbose ); %try Gurobi
            catch;
                error( 'computepolytopenorm:gurobi', 'Gurobi seems to be installed, but not working.' ); end;
        else
            [normval,iteration] = computepolytopenorm_matlab( pts, VV, algorithm, numcore, epsequal, verbose ); end; end; %try matlab
            
end

function [normval,iteration] = computepolytopenorm_matlab( pts, VV, algorithm, numcore, epsequal, verbose );
%test all new points
% pts                   points to test
% VV                    polytope for the norm
% algorithm             type of algorithm
% epsilon               epsilon used in the algorithm
% verbose               verbose level
% norm_lvl              if computed norm is larger than norm_lvl, then tries to compute it exact. 
% removevertexflag      only possible if pts == V. Then, the point pts(:,i) is removed from the polytope prior to computation

    %val=triangulationm(VV);
    %maxdist=max(val);
    %mindist=min(val);

    %point_VV_dist=min(pdist2(pts.',VV.'), pdist2(pts.',-VV.')).';
    %VV_VV_dist= min(pdist2(VV.',VV.'),pdist2(VV.',-VV.'));
    
    %fprintf('Maxdist: %i, ',maxdist); % DEBUG
    %fprintf('Mindist: %i, \n', mindist);
    %fprintf('Polytope size: %i, \n',size(VV,2)); % DEBUG

    num = size( pts, 2 );
    vprintf( ' ', 'imp',[1 verbose] );
    vprintf( [ '\n' repmat('.',1,num) '\n\n'], 'imp',[2 verbose] ); %make waiting-string
    normval=inf*ones(1,num); %infintiy means outside
    iteration=zeros(1,num);
    
    %p=gcp('nocreate'); %do not create parpool explicitely, since this will not work if the parallel toolbox is not installed
    parfor(i=1:num,numcore) 
    %for i=1:num     
            %[~,idx]=min(point_VV_dist(:,i)); %search nearest point
            %idx=VV_VV_dist(:,idx)<1.2*maxdist;  %take all points near to that point
            
            %[normval(i), iterations(i)]=tjsr_normfunct_matlab(pts(:,i), VV(:,idx), 1e-9, verbose-2, algorithm);

             %if(normval(i)>norm_lvl)
                 [normval(i),iter2]=normfunct_matlab(pts(:,i), VV, 1e-9, verbose-2, algorithm);
                 iteration(i)=iteration(i)+iter2;
                 %fprintf('\b!\n'); 
             %end

            if(isnan(normval(i)));          if(verbose>=1); fprintf('\bE\n'); end; normval(i)=inf; 
            elseif(~isfinite(normval(i)));  if(verbose>=1); fprintf('\b8\n'); end; normval(i)=inf; 
            elseif(normval(i)<0);           if(verbose>=1); fprintf('\bm\n'); end; normval(i)=inf; 
            elseif(normval(i)<1-epsequal);  if(verbose>=1); fprintf('\b_\n'); end;
            elseif(normval(i)<1+eps);       if(verbose>=1); fprintf('\b.\n'); end;
            elseif(normval(i)<1+1000*eps);  if(verbose>=1); fprintf('\b,\n'); end;
            elseif(normval(i)<2);           if(verbose>=1); fprintf('\bo\n'); end;
            else;                           if(verbose>=1); fprintf('\bO\n'); end;
            end; end; %prints '-'/number, whether vertex is inside or outside of the polytope    
    vprintf( '\b', 'imp',[1 verbose] );
    vprintf( '\n', 'imp',[2 verbose] ); 
end

function [m, iter]=normfunct_matlab(p, VV, ~, verbose, algorithm)
    % DEBUG; search if point p is contained in vertices V. If so, remove it
     idx = all(bsxfun(@eq, VV, p), 1);
     VV(:,idx)=[];

    if(isequal(p,zeros(size(VV,1),1))); 
        m=0; 
        iter=0; return; end;
    
    n=size(VV,2); %number of vertices of polytpe
    dim=size(VV,1); %dimension
    
    if(isempty(VV)); 
        m=inf; 
        iter=0; return; end;
    
    SCALE=norm(p,1); p=p/SCALE; %scale the problem
    
    switch algorithm;
        case TJSR_MINKFUNCT;
            f=zeros(1,2*n+1);
            f(end)=-1;
            Aeq=[VV zeros(dim,n) -p]; 
            beq=zeros(dim,1);
            A=vertcat([sparse(1,n) sparse(ones(1,n)) 0], [-speye(n) -speye(n) sparse(n,1)], [speye(n) -speye(n) sparse(n,1)], [sparse(n,n) -speye(n) sparse(n,1)]);
            b=vertcat(1, sparse(3*n,1));
            
        case TJSR_CONEFUNCT;           
            f=zeros(1,n+1); f(end)=-1;
            A=spalloc(dim+1+n,n+1,dim*n+n+n+dim);
            A(:,1:n)=[-VV; sparse(ones(1,n)); -speye(n)];
            A(1:dim,end)=p(:); 
            b=[sparse(dim,1); 1; sparse(n,1)];
            sze=size(A);
            Aeq=sparse(1,sze(2));
            beq=0;
            
        case TJSR_COMPLEXFUNCT;
            error('TJSR_COMPLEXFUNCT not implemented yet.');
            
        otherwise;
            error('Unkown value for ''algorithm''.'); end;
    

    %choose algorithm, interior point algorithm does not work well.
    if( verLessThan('matlab','8.1') ); 
        opts = optimset( 'Display','off', 'LargeScale','off', 'Simplex','on' );  %prior 2013 %%XX Tolerance parameter not tested, if this is the right way to set this parameter
    elseif verLessThan('matlab','8.4'); 
        opts = optimset( 'Display','off', 'Algorithm','simplex' ); %prior 2014a %XX Tolerance parameter not tested, if this is the right way to set this parameter
    else; 
        opts = optimoptions( @linprog, 'Display','off', 'Algorithm','dual-simplex' ); end; %after 2014b 
    
    %opts=[]; %for debugging if necessary
    
    %Solve the problem
    [m,~,exitflag,info]=linprog(f,A,b,Aeq,beq,[],[],[],opts); 
    iter=info.iterations;
    if(isempty(m)); 
        m=NaN; 
    else; 
        m=m(end)/SCALE; end;

    %interpret error codes and return value. Output Debug Information
    if(exitflag<=0 || isnan(m))
        if(verbose>=1); fprintf('E'); end;
        if(verbose>=2); fprintf('rror (minkfunct)\n'); end;
        if(verbose>=2); exitflag, end; %#ok<NOPRT>
        m=NaN; return;
    end

    %Text Output
    if(m==0  && rank(VV)==dim && exitflag~=1 ); %dont repair overflows.
        if(verbose>=1); fprintf('U'); end;
        if(verbose>=2); fprintf('nderflow (normfunct)\n'); end;
        m=Inf; return; end;
    if(m<0);
        if(verbose>=1); fprintf('L'); end;
        if(verbose>=2); fprintf('ess than zero (normfunct).\n'); end; end;
    m=1/m;

end

function [normval,iteration]=computepolytopenorm_gurobi(pts, VV, algorithm, numcore, epsequal, verbose);  %epsilon is unused
%test all new points
% pts                   points to test
% VV                    polytope for the norm
% algorithm             type of algorithm
% epsilon               epsilon used in the algorithm
% verbose               verbose level
% norm_lvl              not used

%    if(algorithm~=1); error('Only minkfunct possible for GUROBI'); end; %minkfunct
    n=size(VV,2); %number of vertices of polytpe
    dim=size(VV,1); %dimension
    p = gcp();
    numworker=p.NumWorkers;
    
    %XX Check if there is only one point, 
    
    %search for minimal spanning tree
    [chain, MST]=partitionatepolytope(pts,ceil(log2(numworker)+1));
    numchains=size(chain,2);
    
    %Initialize variables and construct model/params
    num=cellfun(@(x) x.numnodes, MST);
    nearpointidx=cell(1,numchains); 
    for w=1:numchains; 
        [~,nearpointidx{w}] = min(pdist2(MST{w}.Nodes.co,VV.').',[],1); end 
    
    switch algorithm
        case TJSR_MINKFUNCT
            model.obj=zeros(1,2*n+1); model.obj(end)=1; %in matlab terminology: f
            A=vertcat([sparse(1,n) sparse(ones(1,n)) 0], [-speye(n) -speye(n) sparse(n,1)], [speye(n) -speye(n) sparse(n,1)], [sparse(n,n) -speye(n) sparse(n,1)]);
            Aeq_without_p=[VV sparse(dim,n) zeros(dim,1)];  %Aeq=[V sparse(dim,n) -p]; 
            model.A=[A; Aeq_without_p]; %lhs must be sparse
            b=vertcat(1, zeros(3*n,1)); %rhs must be dense
            beq=zeros(dim,1);
            model.rhs=[b; beq]; %in matlab terminology: b
            model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq_without_p,1),1)];

            model.lb = -inf*ones(size(model.A,2),1);
            %model.lb(end) = 0.2; %we are only interested for norm-values of points near to the boundary
            %model.ub =  2*ones(size(model.A,2),1); 
            model.ub = 100*ones(size(model.A,2),1); %if point to be computed is very small, this prevents the solution from becoming unbounded
            %model.ub(end) = 2; %this does not help %we are only interested for norm-values of points near to the boundary 
            model.modelsense = 'max';
            clear A Aeq_without_p b beq;
            
        case TJSR_CONEFUNCT
            model.obj=zeros(1,n+1); model.obj(end)=-1; %in matlab terminology: f
            A=spalloc(dim+1+n,n+1,dim*n+n+n+dim);
            A(:,1:n)=[-VV; sparse(ones(1,n)); -speye(n)];
            %A(1:dim,end)=p(:); %is done in computepolytopenorm_gurobi_worker
            sze=size(A);
            Aeq=zeros(1,sze(2));
            model.A=[A; Aeq]; %lhs must be sparse
            
            b=[zeros(dim,1); 1; zeros(n,1)];
            beq=0;
            model.rhs=[b; beq]; %in matlab terminology: b
            model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
            model.modelsense = 'min';
            
            clear A b beq;
            
        otherwise
            error('not implemented.'); end;
          

    params.outputflag=tif(verbose>=4,1,0);
    %params.FeasibilityTol = max(min(epsilon,1e-2),1e-9); %smaller epsilon is faster
    %params.FeasibilityTol = 1e-9; %smaller epsilon is faster
    params.Method=0; %0=primal simplex (fastest?), 1=dual simplex, 3=barrier (interior point) (slowest?)   
    %params.Presolve=2; %slower?, useless
    %params.DualReductions=0; %useless
    
    vprintf(['\n' repmat('.',1,sum(num)) '\n'],'imp',[3 verbose]); %make waiting-string%make waiting-string
    
    % parfeval
    %%%%%%%%%%%%%%%%%%%%
    [~,idx] = sort( num, 'descend' ); %sort in descending order, such that the ones which take hopefully the longest start at first. We hope that this leads to a 100% CPU load most of the time.
    if( numcore>1 )
        for w = idx;
            %computepolytopenorm_gurobi_worker(n, chain{w}, model, params, MST{w}, nearpointidx{w}, algorithm, epsequal); %DEBUG
            MSTpool(w) = parfeval(p,@computepolytopenorm_gurobi_worker,1,n, chain{w}, model, params, MST{w}, nearpointidx{w}, algorithm, epsequal); end; %#ok<AGROW>
        

        % Collect the results as they become available.
        for w = numchains:-1:1;
            [idx,value] = fetchNext(MSTpool);  % fetchNext blocks until next results are available.
            MSTresults{idx} = value;
            vprintf('%c', (MSTresults{idx}.Nodes.char).','imp',[2 verbose]); end;
        MST = MSTresults;
        clear MSTresults;
    else
        MSTresults=cell(1,numchains);
        for w=1:numchains;
            MSTresults{w}=computepolytopenorm_gurobi_worker(n,chain{w},model,params, MST{w}, nearpointidx{w}, algorithm, epsequal);
            vprintf('\b%c\n', (MSTresults{idx}.Nodes.char).','imp',[2 verbose]); end; end;

    %%%parfor %Same as above but with a parfor loop. This is slower
    %%%%%%%%%%%%%%%%    
%     model=repcell(model,[1 numchains]);
%     params=repcell(params,[1 numchains]);
%     MSTresults=cell(1,numchains); %parfor loop
%     vprintf('\n','imp',[2 verbose]); %make waiting-string%make waiting-string
%     for w=1:numchains    %% for loop
%     %parfor(w=1:numchains,numworker);
%         MSTresults{w}=computepolytopenorm_gurobi_worker(n, chain{w}, model{w}, params{w}, MST{w}, nearpointidx{w}, algorithm, epsequal);
%         vprintf('\b%c\n', (MSTresults{w}.Nodes.char).','imp',[1 verbose]);
%     end
%     MST=MSTresults;
%     clear MSTresults;
%     
    vprintf('\n','imp',[3,verbose]); 
    iteration = sum(cellfun(@(x) sum(x.Nodes.iterations), MST));
    normval = cell2mat(cellfun(@(x) [double(x.Nodes.idx), x.Nodes.norm] , MST,'UniformOutput',0)');
    normval = sortrows(normval);
    normval = normval(:,2).';
    
% %%%%%%%Plot MST%%%%%%%    
%     clf; hold on;
%     plotm(VV,'x');
%     plotm(-VV,'.');
%     h=cell(1,numchains);
%     for w=1:numchains
%         h{w}=plot(MST{w},'k','LineWidth',1,'Layout','layered');
%         highlight(h{w},1);
%         drawnow;
%         h{w}.XData=MST{w}.Nodes.co(:,1);
%         h{w}.YData=MST{w}.Nodes.co(:,2);
%         if(dim>=3); h{w}.ZData=MST{w}.Nodes.co(:,3); end;
%
%         %h{w}.NodeLabel=MST{w}.Nodes.norm;
%         h{w}.NodeLabel=MST{w}.Nodes.iterations;
%      end
%     axis equal
%     drawnow;
    
    


end

function [MST] = computepolytopenorm_gurobi_worker( n, chain, model, params, MST, nearpointidx, algorithm, epsequal ); 

    num = size( chain, 1 );
    dim = size( MST.Nodes.co, 2 );
    normval = zeros( num, 1 );
    iteration = zeros( num, 1 );
    co = MST.Nodes.co;
    
    switch algorithm
        case TJSR_MINKFUNCT; 
            cbasis = int8( zeros(num,dim+3*n+1) );  %vbasis not needed, since it is always zero
            
        case TJSR_CONEFUNCT; 
            sze = size( model.A );
            cbasis = int8( zeros(num,sze(1)) );  %vbasis not needed, since it is always zero
            vbasis = int8( zeros(num,sze(2)) ); end;

    for i = 1:num
        idxnew = chain(i,2); %idx of point to be processed
        idxold = chain(i,1); %idx of parent vertex
        switch algorithm
            case TJSR_MINKFUNCT; 
                model.A(end-dim+1:end,end) = -MST.Nodes.co(idxnew,:);
                %warm start
                if( idxold ); 
                    % % if point is near to a point which has already been computed
                    % % It seems that a basis from a solution from a real vertex is always better than a manually basis from a vertex-point from the polytope.
                    % % Thus we do not need the distances between the points in the chain for the linear-programming part
                     model.cbasis = double( cbasis(idxold,:).' );
                     model.vbasis = zeros( 2*n+1, 1 );
                else
                    % %if no point near to this point has been computed so far, we make a manual basis stemming from a solution of a vertex point of the polytope
                    idx = nearpointidx(idxnew);
                    cbasisVV = zeros( 3*n+1+dim, 1 );
                    cbasisVV(1) = -1; 
                    cbasisVV(end-dim+1:end) = -1;
                    cbasisVV(1+idx+n) = -1; 
                    cbasisVV(2+2*n:idx+2*n) = -1; 
                    cbasisVV(2+2*n+idx:end-dim) = -1;
                    missing = 2*n+1-nnz( cbasisVV );
                    if( missing>0 ); 
                        val = min( idx-1, missing );
                        cbasisVV(2+n:val+2+n-1) = -6; end;
                    missing = 2*n+1-nnz( cbasisVV );
                    if( missing>0 ); 
                         val = min( n-idx, missing );
                         cbasisVV(2+n+idx:val+2+n+idx-1) = -7; end;
                    model.cbasis = cbasisVV;
                    model.vbasis = zeros(2*n+1,1); end;
                
            case TJSR_CONEFUNCT; 
                model.A(1:dim,end) = co(idxnew,:);
                if( idxold ); 
                    % % if point is near to a point which has already been computed
                    % % It seems that a basis from a solution from a real vertex is always better than a manually basis from a vertex-point from the polytope.
                    % % Thus we do not need the distances between the points in the chain for the linear-programming part
                    model.cbasis = double( cbasis(idxold,:).' );
                    model.vbasis = double( vbasis(idxold,:).' ); end; end;

        result = gurobi( model, params );

        iteration(idxnew) = result.itercount;
        if( ~isequal(result.status,'OPTIMAL') );
            normval(idxnew) = NaN;
        else
            cbasis(idxnew,:) = int8( result.cbasis.' );
            if( algorithm==TJSR_CONEFUNCT ); 
                vbasis(idxnew,:) = int8( result.vbasis.' ); end;
            normval(idxnew) = 1/result.x(end); end; end;
    
    % Text Output
    %%%%%%%%%%%%%%%%%%%%5
    strval = 79*ones( size(normval) );              % O
    %idx = isnan( normval );     strval(idx)=78;     % N %slow
    idx = ~isfinite( normval ); 
    strval(idx) = 56;     % 8
    idx = normval<2;          
    strval(idx) = 111;    % o
    idx = normval<1+1000*eps;      
    strval(idx) = 44;     % ,
    idx = normval<1+eps;      
    strval(idx) = 46;     % .
    idx = normval<1-epsequal;      
    strval(idx) = 95;     % _
    idx = normval<0;          
    strval(idx) = 109;    % m
    
    % Save Output
    %%%%%%%%%%%%%%%5
    normval( ~isfinite(normval) ) = inf;
    normval( normval<0 ) = inf;
    
    MST.Nodes.norm = normval;
    MST.Nodes.char = strval;
    MST.Nodes.iterations = iteration;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 