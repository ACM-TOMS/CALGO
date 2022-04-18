function [pn_est, region] = estimatepolytopenorm(pts, VV, VV2, ptype, epsilon);
% [pn_est] = estimatepolytopenorm(pts, VV, VV2, [ptype], epsilon);
% Estimates the polytope_norm of pts wrt VV
%
% Input:
%   pts             points
%   VV              vertices of the polytope. 
%   VV2             vertices of the polytope and vertices which are inside of the polytope. Can be empty.
%   ptype           scalar, defines the used estimate/type of polytope w.r.t to given vertices. It is one of the following values
%                       0: positive cone (for cone-norm)
%                       1: (default) point symmetric (for minkowski-norm), vertices of the polytope are [VV -VV];
%                       2: complex polytope (for complex-polytope-norm)
%                       otherwise: undefined behaviour
%   epsilon         Epsilon used to determine whether a point is inside or outside
%               
% Output:
%   pn_est          polytope-norm estimate. Computation is based on type:
%   region          1x2 array. Gives intervals where the point is for sure inside/outside the polytope
%                       For region = [a b]
%                       if 0 <= pn_est <  a : p is inside
%                       if a <= pn_est <= b : nothing is known
%                       if b <  pn_est <= inf : p is outside
%
% Notes:
%        If pts is empty pn_est is empty
%        If VV is empty, pn_est=inf
%        This function is optimized to handle a lot of points p at once
%        There are more ptypes available, but they are not tested at all
%                       5: experimental minkowski-norm. Searches for the midpoints of the nearest faces. works only for dim 2
%                       6: diff-distance wrt to [VV -VV]
%                       10: conefunct
%                       11: minkfunct
%
% Notes for ptype==0:
%       all entries in p and VV must be nonnegative
%       If pn_est<1, then the norm of the point is proofed to be inside the cone!
%       If the set of points in VV is dense, then the estimate is also not bad for points outside
%
% E.g.: normval=estimatepolytopenorm(randn(3,10), randn(3,10), [], 11);
%
% See also: testnorm
%
% Written by tommsch, 2018

% Changelog: 2019_10_08 - tommsch - Code for type 0 (positive cone) rewritten
%                                   Added option "epsilon"

% XX Return true bound for norm
% XX Use epsilon to determine inside/outside. Only done for type 0 yet.
% XX Implement search for nearest face also for case (P)
% XX think about estimate_outside_4 for case (P)

if(nargin~=5); fprintf('Wrong number of input arguments'); error('Error'); end;



if(isempty(pts)); 
    region=[0 inf];
    pn_est=[]; 
    return;
elseif(isempty(VV)); 
    region=[0 inf];
    pn_est=inf(1,size(pts,2)); 
    return;
end


switch ptype
    case {0,'defaultP','defP'}; %TJSR_CONEFUNCT
        
        %norm values for inside/outside
        IN=.5;
        OUT=2;
        
        sze=size(VV,2);
        sze2=size(VV2,2)+size(VV,2);
        nump=size(pts,2);
        [estimate_inside]=deal(zeros(1,nump));
        [estimate_outside1,estimate_outside2,estimate_outside3]=deal(ones(1,nump));
        %estimate_outside4=ones(1,nump);
        normpt=sum(pts.^2,1);
        for i=1:nump;
            pt=pts(:,i);
            %Count in how many hypercubes defined by VV the point is not. If at least outside one, then the point is inside.
            estimate_inside(i) = nnz(any([VV VV2]<=pt+epsilon,1));
            %if(estimate_inside(i)==sze2); %if point is not inside, we check if it is outside
                %Project all vertexpoints onto the vector defined by p. If all Vertex points lie on one side of p, then p is outside
                estimate_outside1(i) = nnz(pt'*VV < normpt(i)-epsilon);
                %if one coordinate of the point is larger than all corresponding coordinates of the polytope, then p is outside
                estimate_outside2(i) = nnz(all(pt-epsilon>VV,2));
                %distance to nearest point
                estimate_outside3(i) = min(sum(abs(VV-pt)));
                %count in how many anti-cubes the point lies
                %this is only a true characterization if there are no interior points in VV
                %estimate_outside4(i) = nnz(any([VV VV2]>=pt-epsilon,1));
            %end
        end
        
        %indices of point outside, inside and unkown
        %idx_outside = estimate_outside1==sze | estimate_outside2>0 | estimate_outside4~=sze2;
        idx_outside = estimate_outside1==sze | estimate_outside2>0;
        idx_inside = estimate_inside~=sze2;
        idx_unkown = ~(idx_outside | idx_inside);
        
        %normalize values
        estimate_outside1=estimate_outside1./sze;
        estimate_inside=(sze2-estimate_inside)./sze2;
        %estimate_outside4=1+(sze2-estimate_outside4)./sze2;
        
        %set norms
        val=1-estimate_inside;
        pn_est(idx_inside)=min(IN-.1,val(idx_inside));;
        val=1/2+2*estimate_outside3+2*estimate_outside1;
        pn_est(idx_outside)=max(OUT+.1,val(idx_outside));  
        val=1/2+2*estimate_outside3+2*estimate_outside1-estimate_inside;
        pn_est(idx_unkown)=max(IN+.1,min(val(idx_unkown),OUT-.1));        

        region=[IN OUT];            

    case {1,'defaultR','defR'} %TJSR_MINKFUNCT
        
        LOW=1;
        HIGH=2;
        sze=size(VV,2);
        dim=size(VV,1);
        
        %norm of the point corresponding to an inscribed rectangle
         if(rank(VV)>=dim)
            VVrand = randpart(VV,'size',5*dim,'minsize',2*dim,'maxnum',10*dim,'cols');
            for i=1:size(VVrand,2)
                if(rank(VVrand{i})<dim); VVrand{i}=[]; end; %remove sets which have not full rank
            end
            VVrand(cellfun(@isempty,VVrand))=[];
            VVrand{end+1}=VV; %to make sure that there is at least on thing there

            val=inf;
            try %here we have memory problems if VVrand is really big
                for i=1:size(VVrand,2);
                    val=min(val, sum(abs(VVrand{i}\pts),1));
                end
            catch
                fprintf('''estimatepolytopenorm'': Error (Probably out of memory).\n');
            end
                
            estimate_inside=val;

         else
            estimate_inside=inf*ones(1,size(pts,2));
         end
        
        
        idx=estimate_inside<1; %it is important that here is a "less". Otherwise, points which are equal to vertex-points are reported as inside
        pts(:,idx)=[]; %test only points which are not inside, if they are outside
        try %against OUT OF MEMORY errors
            val=pts'*[VV -VV];
            normpt=sum(pts.^2,1);
            estimate_outside1=val<normpt';
            estimate_outside1=sum(estimate_outside1,2);
        catch
            estimate_outside1=zeros(size(pts,2),1);
        end
        
        
        pn_est=min(HIGH,estimate_inside); 
        
        idx2=estimate_outside1==sze*2;
        estimate_outside1(idx2)=estimate_outside1(idx2)/(sze*2)*(HIGH+0.1);
        estimate_outside1(~idx2)=max(LOW,estimate_outside1(~idx2)/(sze*2)*(HIGH-0.1));
        pn_est(~idx)=estimate_outside1;
        
        region=[1 2];
        

        
        
    case {5,'project2'}; %TJSR_MINKFUNCT experimental %works only for 2d
        %searches a good vector to project on
        %VV=[VV -VV];
        R = 2; 
        N = 500;
        spherepts = linspace(0,2*pi,N+1); 
        spherepts = R*[cos(spherepts); sin(spherepts)]; 
        spherepts(:,end)=[];
        nspherepts=size(spherepts,2);
        nVV=size(VV,2);
        dim=size(VV,1);
        if(~isequal(dim,2));
            error('ptype 5 works only for dimension 2. Use ptype 55 instead.\n'); end;
        
        npts=size(pts,2);
        pn_est=zeros(1,npts);
        
        for k=1:npts
            
            
%             %weighted distances
%             POWER=100;
%             weightspherepts=zeros(1,nspherepts);
%             for i=1:nspherepts
%                 sphereptsi=spherepts(:,i);
%                 weightspherepts(i)=norm(sphereptsi-pts(:,k))^POWER*nVV;
%                 weightspherepts(i)=weightspherepts(i)-sum(sum((VV-sphereptsi).^2,1).^(POWER/2));
%             end

            [~,idx]=min(sum((spherepts-pts(:,k)).^2,1).^(1/2)); %search nearest spherept
            if(idx>=nspherepts);
                val1=diffdistance_worker(spherepts(:,1),VV);
            else
                val1=diffdistance_worker(spherepts(:,idx+1),VV);
            end
            val0=diffdistance_worker(spherepts(:,idx),VV);
            if(val1<val0); 
                dir=+1; 
            else; 
                dir=-1; 
                %val0=val1; 
            end; % i am looking for a local minimum
            while(true)
                idx=idx+dir;
                if(idx==0); idx=nspherepts; end;
                if(idx>nspherepts); idx=1; end;
                val1=diffdistance_worker(spherepts(:,idx),VV);
                if(val1<val0); 
                    val0=val1;
                else
                    idx=idx-dir;
                    break;
                end;
            end
            if(idx==0); idx=nspherepts; end;
            if(idx>nspherepts); idx=1; end;
            dir=spherepts(:,idx)/R;
            
            pn_est(k)=nnz(dir'*[VV -VV]<dir'*pts(:,k));
        end
        idx=pn_est==2*nVV;
        pn_est(idx)=3;
        pn_est(~idx)=pn_est(~idx)/(2*nVV);
        region=[0 2];
        
    case {55,'project'}; %TJSR_MINKFUNCT experimental
        %searches a good vector to project on,
        %should work for any dimension
        %VV=[VV -VV];
        R=5;
        nVV=size(VV,2);
        
        npts=size(pts,2);
        pn_est=zeros(1,npts);
        
        parfor k=1:npts
            x0=pts(:,k);
            SCOx0=cart2sphm2(x0);
            fun = @(x) diffdistance_worker( sph2cartm2(x)'*R, VV );
            
            SCOx = fminsearch(fun,SCOx0);
            if(isnan(SCOx));
                2;
                SCOx=SCOx0;
            end
            dir=sph2cartm2(SCOx)';
            
            pn_est(k)=nnz(dir'*[VV -VV]<dir'*pts(:,k));
        end
        idx=pn_est==2*nVV;
        pn_est(idx)=3;
        pn_est(~idx)=pn_est(~idx)/(2*nVV);
        region=[0 2];
             
        
    case {6,'diffdistance','diffdist'} %diffdistance
        npts=size(pts,2);
        nVV=size(VV,2);
        pn_est=zeros(1,npts);
        VV=VV.';
        NORM=sum(VV.^2,2).^(1/2);
         for k=1:npts
             p=pts(:,k).';
              d1=pdist2(p,VV).'; %transpose since matlab stupidly returns row-vector
              d2=pdist2(p,-VV).';
              %d1=((d1)./NORM);
              %d2=((d2)./NORM);
              pn_est(k)=max(abs(d1-d2));
             
%             prep=repmat(p,[nVV 1]);
%             pn_est(k)=max(abs(sum(prep.*VV,2)));
             
             
         end
         
         region = [0 inf];
         
    case {7,'ellipse','ell'} %ellipse
        
        l=pinv(VV)*pts;
        p=2; pn_out=sum(abs(l).^p,1).^(1/p); %p must be 2
        
        p=1; pn_in1=sum(abs(l).^p,1).^(1/p); %p must be 1
        l=VV\pts;
        p=1; pn_in2=sum(abs(l).^p,1).^(1/p); %p must be 1
        
        pn_in=min(pn_in1,pn_in2);
        
        idx_in=pn_in<1;
        idx=idx_in;
        pn_est(idx)=min(1*pn_in);
        
        idx_out=pn_out>1;
        idx=idx_out;
        pn_est(idx)=2*pn_out(idx);
        
        idx=~idx_out & ~idx_in;
        pn_est(idx)=max(1,min(2,(pn_out(idx) +1*pn_in(idx))/2));
        
        region = [ 1 2 ];        

    case {2,'defaultC'}; %TJSR_COMPLEXFUNCT;
        error('not implemented yet.');
        
    case {10, 'P'};
        pn_est = computepolytopenorm( pts, VV, 0, 0, eps, 0 );
        region = [1 1];
        
    case {11, 'R'};
        pn_est = computepolytopenorm( pts, VV, 1, 0, eps, 0 );
        region = [1 1];
        
    case {12, 'C'};
        error('not implemented yet.');
        
    otherwise;
        error('Wrong ptype given.');

end;

end

function dm = diffdistance_worker(p, VV);
    p=p.'; VV=VV.';
    d1=pdist2(p,VV).'; %transpose since matlab stupidly returns row-vector
    d2=pdist2(p,-VV).';
    NORM=sum(VV.^2,2).^(1/2);

    %misst Abstand 
    %d=sort(min(d1,d2));
    %if(length(d)>=1); d=d(1); else; d=inf; end;


    %findet Punkte die weit weg von anderen Punkten sind UND die innerhalb sind
    d1=(d1./NORM);
    d2=(d2./NORM);
    dm=max(abs(d1-d2));
    
        

end

function ret = randpart(varargin);
% [k | s] = randPartition(K | S, ['size',val | 'num',val], [options]);
% divides a finite set of linear size K into n, equally big, partitions.
%
% Input:
%   K               (integer scalar) number of elements in the set, or
%   S               (not a scalar) the set which we want to partitionate
%   'size',val      returns equally sized partitions of roughly size val
%   'num',val       returns num many partitions of roughly equal size
%                   Undefined behaviour for val not being an integer between 1 and the number of elements
%
% Options:
%   'maxnum',val    Maximum number of partitions
%   'minsize',val   Minimum size of partitions
%                   Undefined behaviour if 'size','num','maxnum' and 'minsize' are not compatible
%   'cols'          Partitionates the columns of the matrix S.
%   'rows'          Partitionates the rows of the matrix S.
%   'verbose',val   Verbose level
%   
% Output:
%   k       cell array of linear index arrays OR
%   s       cell array of partitions
%
% E.g.: 
%
% See also: randperm
%
% Written by tommsch, 2018

if(isscalar(varargin{1}) && iswholenumber(varargin{1})); 
    nS=varargin{1}; %number of elements in S
else; 
    S=varargin{1};
    if(parsem('cols',varargin));
        nS=size(S,2);
    elseif(parsem('rows',varargin));
        nS=size(S,1);
    else
        nS=numel(varargin{1});
    end
end

maxnum=parsem('maxnum',varargin,inf);
minsize=parsem('minsize',varargin,0);
if(minsize>nS); minsize=nS; end;


if(parsem('num',varargin))
    num=varargin{3}; %number of partitions
else %case: 'size'
    num=round(nS/varargin{3});
end
if(num < 0 || ~isfinite(num)); error('Value for ''size'' or ''num'' is bad.'); end;
if(num== 0); num=1; end;
if(num>=nS); num=nS; end;
if(num>maxnum); num = maxnum; end;
if(minsize>nS/num); 
    num=ceil(nS/minsize); 
end;
while(true)
    edges = round(linspace(1,nS+1,num+1));
    if(min(diff(edges))>=minsize); break; else; num=num-1; end;
    if(num<=0); error('Value for ''minsize'' or ''maxnum'' is bad.'); end;
end

k=size(edges,2)-1; %number of partitions

pos = randperm(nS);

ret = cell(1,k);
for i = 1:k
    ret{i} = edges(i):edges(i+1)-1;
    ret{i}=pos(ret{i});
end

if(isscalar(varargin{1}) && iswholenumber(varargin{1})); 
    %do nothing
else
    if(parsem('cols',varargin));
        for i = 1:k; ret{i}=S(:,ret{i}); end;
    elseif(parsem('rows',varargin));
        for i = 1:k; ret{i}=S(ret{i},:); end;
    else
        for i = 1:k; ret{i}=S(ret{i}); end;
    end
end


end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   
