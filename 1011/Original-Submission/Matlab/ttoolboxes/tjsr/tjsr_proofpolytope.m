function type=tjsr_proofpolytope(M, type, prooftype); 
% type = tjsr_proofpolytope(M, type, prooftype); 
% This function belongs to tjsr!
%       This function uses (nearly) no tricks to accelerate the computation
%       Every trick is marked with the label %TRICK
%       It uses matlabs (slow) linprog to proof the polytope, since here mistakes are more unlikely
%       It stores its output in the struct type.proof
%       If prooftype==2, then the normalized/preprocessed matrices are used.
%
% Input:
%       M               the original (unscaled,...) matrices
%       type            the output from tjsr
%       prooftype       either 1 or 2. Determines which matrices are used to proof the polytope
%
% Output:
%       type.proof.V            the polytope used in the proof
%       type.proof.M            the matrices used in the proof
%       type.proof.pts          the points used in the proof (pts=[M{i}*V])
%       type.proof.algorithm    the algorithm used
%       type.proof.norm         the computed norms. 
%       type.proof.errorflag    
%                               -2   if there is nothing to proof
%                               -1   if the proof terminated without problems and the polytope is likely to be invariant
%                               0    if the proof terminated without problems, but the polytope is likely NOT to be invariant
%                               1    if the proof terminated early
%                               2    if we could not proof for some specific reasons
%                               nan  if something strange happened
%       type.proof.JSR          JSR estimate from proof
%
% Written by tommsch, 2018

    type.proof.errorflag=nan;
    try    
        if(~isfield(type,'cyclictree')); 
            type.proof.errorflag=-2;
            type.info.infotext=vprintf('Cannot proof polytope, since there is no field ''cyclictree'' inside the struct ''type.''\n. Either the dimension of the matrices is 1, there are invariant subspaces or something else happened.\n','imp',[1,type.opt.verbose],'str',type.info.infotext);
            return;
        elseif(isfield(type,'block'));
            type.proof.errorflag=2;;
            type.info.infotext=vprintf('Cannot proof polytope, since the matrices are block-diagonalized.\n','imp',[1,type.opt.verbose],'str',type.info.infotext);
            return;
        else
            type.info.infotext=vprintf('Proof invariance of polytope. Time: %s\n',datetime('now'),'imp',[1,type.opt.verbose],'str',type.info.infotext);
        end


        VV = [type.cyclictree.V{:}];
        type.proof.VV=VV;

        if(prooftype==1);
            JSR=type.lambda;
            %do nothing
        else
            JSR=1;
            M=type.M_normalized;
        end
        type.proof.M=M;

        pts=cell(1,type.counter.nummatrix);
        for i=1:type.counter.nummatrix
            pts{i}=M{i}*VV;
        end
        pts=cell2mat(pts);
        type.proof.pts=pts;

        if(type.info.algorithm==TJSR_CONEFUNCT)
            type.proof.algorithm=TJSR_CONEFUNCT;
            normval=conefunct(pts,VV,type.opt.verbose,JSR);
        elseif(type.info.algorithm==TJSR_MINKFUNCT)
            type.proof.algorithm=TJSR_MINKFUNCT;
            normval=minkfunct(pts,VV,type.opt.verbose,JSR);
        elseif(type.info.algorithm==TJSR_COMLPEXFUNCT)
            error('not implemented yet');
        end
        type.proof.norm=normval;

        if(prooftype==1)
            MAX=max(normval,[],'includenan');
            m=MAX/type.lambda-1;
        else
            MAX=type.lambda*max(normval);
            m=max(normval)-1;
        end


        type.info.infotext=vprintf('Estimate for the JSR from Proof: [ %15.12g, %15.12g]\n',type.lambda, MAX,'imp',[1,type.opt.verbose],'str',type.info.infotext);
        if(m<1e-7);   
            type.proof.errorflag=-1;
            type.info.infotext=vprintf('Polytope may be invariant. m=1+%e\n',m,'cpr',[0,.5,0],'imp',[1,type.opt.verbose],'str',type.info.infotext);
        else; 
            type.proof.errorflag=0;
            type.info.infotext=vprintf('Polytope is not invariant with m=1+%e\n',m,'cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext);
        end
        type.proof.JSR=[type.lambda MAX];
        
    catch
         type.proof.errorflag=1;
         type.info.infotext=vprintf('Proof terminated early.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext);
    end
end

function [normval]=conefunct(pts,VV,verbose,JSR)
%computes the norm needed for prot2012, Algorithm (P)

    dim=size(VV,1);
    npts=size(pts,2);
    nV=size(VV,2);
    
    if(min(max(VV,[],1))<10*eps)
        vprintf('Warning: Polytope degenerated.\n','cpr','err','imp',[2 verbose]);
    end

    f=zeros(1,nV+1);
    f(end)=-1;
    A=spalloc(dim+1+nV, nV+1, dim*nV+nV+nV+dim); %last column contains p:
    A(:,1:nV)=[-VV; sparse(ones(1,nV)); -speye(nV)];
    b=[sparse(dim,1); 1; sparse(nV,1)];
    try
        if verLessThan('matlab','8.1'); opts=optimset('Display','off','LargeScale','off', 'Simplex','on');  %prior 2013
        elseif verLessThan('matlab','8.4'); opts=optimset('Display','off','Algorithm','simplex'); %prior 2014a
        else; opts = optimoptions(@linprog,'Display','off','Algorithm','dual-simplex'); %after 2014b
        end
        %opts=[]; %for debugging if necessary
    catch
        opts=[];
    end;
    
    %Solve the problem
    normval=zeros(1,npts);
    if(verbose>=1);
        fprintf(['\n' repmat('.',1,npts) '\n\n']);
    end

    parfor i=1:npts
        p=pts(:,i);

        if(isequal(p,zeros(dim,1)));  %
            normval(i)=0; 
        else
            try
                AA=A;
                AA(1:dim,end)=p; %#ok<SPRIX>
                if(isempty(opts))
                    [m,~,exitflag,~]=linprog(f,AA,b);
                else
                    [m,~,exitflag,~]=linprog(f,AA,b,[],[],[],[],[],opts);
                end
                if(exitflag==1);
                    normval(i)=1/m(end);
                elseif(exitflag==-3);
                    normval(i)=inf;
                else
                    normval(i)=nan;
                end
            catch
                normval(i)=nan;
            end
        end
            
        if(verbose>=1)
            val=normval(i)/JSR;
            if(val==inf); fprintf('\b8\n');
            elseif(val>2); fprintf('\bO\n');
            elseif(val>1+100*eps); fprintf('\b;\n');
            elseif(val>1+10*eps); fprintf('\b,\n');
            elseif(val>1-eps); fprintf('\b.\n');
            elseif(val>0); fprintf('\b_\n');
            elseif(val>-inf); fprintf('\bm\n');
            else; fprintf('\b?\n');
            end;
        end

    end

end


function normval=minkfunct(pts, VV, verbose, JSR)
% normval=minkfunct(pts, V, verbose)
% Computes the Minkowski Norm corresponding to polytopes convex hull of the points given by [V -V];
% Thus minkfunct(p,V) and minkfunct(p,[V -V]) give the same number
% The function prints a warning if the computed norm is negative or if the norm could not be computed.
% Each column in V is one vertex of the polytope

    dim=size(VV,1); %dimension
    npts=size(pts,2);
    nV=size(VV,2); %number of vertices of polytpe
    
    if(rank([VV -VV])<dim); 
        vprintf('Warning: Polytope degenerated.\n','cpr','err','imp',[2 verbose]);
    end;

    

    f=zeros(1,2*nV+1);
    f(end)=-1;
    beq=zeros(dim,1);
    A=vertcat([sparse(1,nV) sparse(ones(1,nV)) 0], [-speye(nV) -speye(nV) sparse(nV,1)], [speye(nV) -speye(nV) sparse(nV,1)], [sparse(nV,nV) -speye(nV) sparse(nV,1)]);
    b=vertcat(1, sparse(3*nV,1));
    Aeq=[VV zeros(dim,nV)]; %in the last column we have to add -pts : Aeq=[V zeros(dim,nV) -p];    

    if verLessThan('matlab','8.1'); opts=optimset('Display','off','LargeScale','off', 'Simplex','on');  %prior 2013
    elseif verLessThan('matlab','8.4'); opts=optimset('Display','off','Algorithm','simplex'); %prior 2014a 
    else; opts = optimoptions(@linprog,'Display','off','Algorithm','dual-simplex'); %after 2014b 
    end
    %opts=[]; %for debugging if necessary


    %Solve the problem
    normval=zeros(1,npts);
    if(verbose>=1);
        fprintf(['\n' repmat('.',1,npts) '\n\n']);
    end

    parfor i=1:npts
        p=pts(:,i);    

        if(isequal(p,zeros(dim,1)));  %test if point is zero. Then, the norm is zero.
            normval(i)=0; 
        else
            try
                [m,~,exitflag,~]=linprog(f,A,b,[Aeq -p],beq,[],[],[],opts); 
                if(exitflag==1); 
                    normval(i)=1/m(end);
                elseif(exitflag==-3);
                    normval(i)=inf;
                else
                    normval(i)=nan;
                end  
            catch

                normval(i)=NaN; 
            end
        end

        if(verbose>=1)
            val=normval(i)/JSR;
            if(normval(i)==inf); fprintf('\b8\n');
            elseif(val>2); fprintf('\bO\n');
            elseif(val>1+100*eps); fprintf('\b;\n');
            elseif(val>1+10*eps); fprintf('\b,\n');
            elseif(val>1-eps); fprintf('\b.\n');
            elseif(val>0); fprintf('\b_\n');
            elseif(val>-inf); fprintf('\bm\n');
            else; fprintf('\b?\n');
            end;
        end

    end

end



function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 