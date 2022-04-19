function [triang, blocks, B] = jointTriangul(M,varargin)
%
% JOINTTRIANGUL  Heuristic for joint block-triangularization
%                                        
%
% [TRIANG, BLOCKS, B] = JOINTTRIANGUL(M)
%   If the set M can be jointly block-triangularized then this function 
%   tries to find and return the sets of diagonal blocks. For this it tries     
%   to generate an invariant subspace by applying the matrices in M to each
%   eigenvector of every matrix of M.   
%   
% [TRIANG, BLOCKS, B] = JOINTTRIANGUL(M,V0)
%   Does the same but trying to expand the subspace spanned by the linearly 
%   independent columns of V0 instead of the eigenvectors of the
%   matrices in M.
%
% This is a heuristic and gives no guarantee of finding 
% every jointly invariant subspace.
%
%  TRIANG is a boolean indicating if the matrices could be jointly 
%  triangularized. 
%
%  B is unitary and such that for i=1,...,length(m), if q blocks were
%  found, 
%                  
%                  blocks{1}{i}    *          *       *    . . .    * 
%                     0    
%                              blocks{2}{i}   *       *             .
%                     0            0                                .
%                     .            .       blocks{3}{i}             .
%      B'*M{i}*B =    .            .          0         .           
%                     .            .                        .  
%                                                              .    *
%                     0            0          0     . . .       blocks{q}{i}  
%
%
%  When this is possible, the joint spectral radius of M is equal to the
%  max of the joint spectral radii of each set of BLOCKS.
%  
%
% REFERENCES
%    R.Jungers, 
%      "The Joint Spectral Radius: Theory and Applications" 
%      Vol. 385 section 1.2.2.5 in Lecture Notes in Control and Information
%      Sciences, Springer-Verlag. Berlin Heidelberg, June 2009


% THIS ALGORITHM COULD BE OPTIMIZED
% By putting a flag on the subspaces that have no chance of being 
% re-divided for instance.
% (Tree view)

n = size(M{1},1);
m = length(M);

if n==1
    triang=0;
    blocks=M;
    B = 1;
    
    return;
end

if (nargin==2)
    V0 = varargin{1};    
    [triang, blocks, B] = subJointTriangul(M,V0);    
else    
    triang = 0;
    blocks = {M};
    B = eye(n);
    imat=0;
    
    while (~triang && imat<m)
        imat=imat+1;
        if (issparse(M{imat}) && n>50)
            [V, D] = eigs(M{imat},n);
        else
            [V, D] = eig(full(M{imat}));
        end

        ieigvec = 0;
        while (~triang && ieigvec<n)
            ieigvec=ieigvec+1;
            
            [triang, blo, Bas] = subJointTriangul(M,V(:,ieigvec));
            %fprintf('\nComputing for eigvec %d of matrix %d',ieigvec,imat) 
            
            if (triang)
                blocks = blo;
                B = Bas;
            end

        end

    end
end
end

function [triang, blocks, B] = subJointTriangul(M,V0,step,indnew)
% subJointTriangul(M,V0)
%       
% subJointTriangul(M,V0,step,indnew)
% 
n = size(M{1},1);
m = length(M);

if n==1
    triang=0;
    blocks={M};
    B = 1;
    
    return;
end

if (nargin<3)
    step = 1;
    indnew = 1;
end

B = zeros(n,n);

nv = rank(V0); % Number of linearly independent vectors currently found

[nvnew, Vvecnew] = appMat(M,V0,nv,indnew);

if (nvnew==0 && cond(V0)<1e3)
        
  [Qr,R] = qr(V0);
  Q1 = Qr(:,1:nv);
  Q2 = Qr(:,nv+1:n);
  
  Mtri1 = cell(1,m);
  Mtri2 = cell(1,m);
  
  for imat=1:m
      Mtri1{imat} = Q1'*M{imat}*Q1;
      Mtri2{imat} = Q2'*M{imat}*Q2;
  end
  
  [sub_triang1, sub_blocs1, sub_B1] = jointTriangul(Mtri1);
  [sub_triang2, sub_blocs2, sub_B2] = jointTriangul(Mtri2);
  
  if sub_triang1
      blocks = sub_blocs1;
      B(:,1:nv) = Q1*sub_B1;
  else
      blocks = {Mtri1};
      B(:,1:nv) = Q1;
  end
  
  if sub_triang2
      blocks = [blocks, sub_blocs2];
      B(:,nv+1:n) = Q2*sub_B2;
  else 
      blocks = [blocks, {Mtri2}];
      B(:,nv+1:n) = Q2;
  end
  
  Q.Q1 = Q1;
  Q.Q2 = Q2;
  triang=1;

elseif (nv+nvnew == n || cond(V0)>=1e3)
    if step==1
        V0 = zeros(n,n);
        V0(:,1) = Vvecnew(:,end);
        [triang, blocks, B] = subJointTriangul(M,V0,2,1);
    else % Step == 2, we have generated the whole space from both ends
        triang = 0;
        blocks = M;
        B = eye(n);  
    end
else
    V0 = Vvecnew;
    [triang, blocks, B] = subJointTriangul(M,V0,step,nv+1);
end
end


function [nvnew, Vvecnew] = appMat(M,V0,nv,indNew)
% 
%
%

m = length(M);

Vvecnew = V0(:,1:nv);
nvnew = 0;

for imat = 1:m
    for ivec = indNew:nv
        vectemp = M{imat}*V0(:,ivec);
        normvec = norm(vectemp);
        
        if(normvec > 1e-14)
            Vvectemp = [Vvecnew(:,1:(nv+nvnew)) vectemp/normvec];
            rankTemp = rank(Vvectemp,1e-12);
            
            if (rankTemp == nv+nvnew+1 )
                Vvecnew = Vvectemp;
                nvnew = nvnew+1;
            end
        end
    end
end

end



