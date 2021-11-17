function Qs = mg_smooth_step(As,level,sweeps,smooth,stype)
%mg_smooth_step   smoothers for GMG on step domain
%   Qs = mg_smooth_step(As,level,sweeps,smooth,stype)
%   input
%           As        structure containing coefficient matrices
%           level     finest level of grid
%           sweeps    number of directions used for Gauss-Seidel
%           smooth    choice of smoother (1/Jacobi, 2/Gauss-Seidel, 3/ILU)
%           stype     type of smoother (1/point, 2/line)
%   output
%           Qs        structure containing smoothing operator in factored
%                     form
%
%   IFISS function: AR, HCE; 15 April 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage
for i=level:-1:2;
   A = As(i).matrix;
   N=size(A,1);

   % permute A to lexicographic ordering to help find smoother
   np=2^i; np2=np/2;
   nx1=np2; nx2=5*np2+1; ny1=np2+1;
   nb1=nx1*ny1; norigin=nx2*np2+nb1;
   index=nb1+[1:nx2*np2];
   for row=1:ny1
      col=1:nx1;
      index=[index (row-1)*nx1+col];
      col=1:nx2;
      index=[index norigin+(row-1)*nx2+col];
   end
   nnod=(11/4)*np*np+4*np+1;
   Pmat=sparse(nnod,nnod);
   for ip=1:nnod;
      Pmat(ip,index(ip))=1.0;
   end
   A = Pmat*A*Pmat';

   % number of elements on each dimension of full domain
   nel = 2/3*(sqrt(1+3*N)-2);
   if stype==2
      % line Gauss-Seidel
      % parameters to identify non-zero diagonals
      d1=nel/2+1;d2=nel+1;
      % horizontal lines, sweep bottom to top
      Q1 = tril(A,1);
      [L1,U1] = lu(Q1);
      if sweeps>=2
         % vertical lines, sweep left to right
         Q2 = diag(diag(A,0))+diag( diag(A,-1),-1)+  ...
              diag(diag(A,-d1-1),-d1-1) + diag(diag(A,-d1),-d1)+ ...
              diag(diag(A,-d2-1),-d2-1) + diag(diag(A,-d2),-d2)+ ...
              diag(diag(A,d1-1),d1-1) + diag(diag(A,d1),d1)+ ...
              diag(diag(A,d2-1),d2-1) + diag(diag(A,d2),d2);
         [L2,U2] = lu(Q2);
         if sweeps>=3
            % horizontal lines, sweep top to bottom
            Q3 = triu(A,-1);
            [L3,U3] = lu(Q3);
            if sweeps==4
               % vertical lines, sweep right to left
               Q4 = diag(diag(A,0))+diag( diag(A,1),1)+  ...
                    diag(diag(A,-d1+1),-d1+1) + diag(diag(A,-d1),-d1)+ ...
                    diag(diag(A,-d2+1),-d2+1) + diag(diag(A,-d2),-d2)+ ...
                    diag(diag(A,d1+1),d1+1) + diag(diag(A,d1),d1)+ ...
                    diag(diag(A,d2+1),d2+1) + diag(diag(A,d2),d2);
                    [L4,U4] = lu(Q4);
            else
               L4=sparse(N,N); U4=L4;
            end
         else
            L3=sparse(N,N); U3=L3; L4=L3; U4=L3;
         end
      else
         L2=sparse(N,N); U2=L2; L3=L2; U3=L2; L4=L2; U4=L2;
      end
   else
      % point smoothers
      if smooth==3
         % ILU
         [L1,U1]=ilu0(A);
         L2=sparse(N,N); U2=L2; L3=L2; U3=L2; L4=L2; U4=L2;
      elseif smooth==2
         % point Gauss-Seidel
         Q1 = tril(A,0); [L1,U1] = lu(Q1);
         L2=sparse(N,N); U2=L2; L3=L2; U3=L2; L4=L2; U4=L2;
      else
         % point Jacobi
         omega=8/9; % relaxation factor for damped Jacobi
         Q1 = (1/omega)*spdiags(diag(A),0,N,N);[L1,U1] = lu(Q1);
         L2=sparse(N,N); U2=L2; L3=L2; U3=L2; L4=L2; U4=L2;
      end
   end
%
   % re-permute factors to match original numbering
   L1=Pmat'*L1*Pmat; U1=Pmat'*U1*Pmat; 
   L2=Pmat'*L2*Pmat; U2=Pmat'*U2*Pmat; 
   L3=Pmat'*L3*Pmat; U3=Pmat'*U3*Pmat; 
   L4=Pmat'*L4*Pmat; U4=Pmat'*U4*Pmat; 
%
   Qs(i) = struct('L1',L1,'L2',L2,'L3',L3,'L4',L4, ...
                  'U1',U1,'U2',U2,'U3',U3,'U4',U4);
end
