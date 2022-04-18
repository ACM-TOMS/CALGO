function M = waveletMat(N)
%
%  waveletMat   constructs the set of two matrices whose joint spectral 
%               radius is closely related to the H�lder exponent of 
%               continuity of Daubechies' wavelets [1], [2], [3].
%
% M = waveletMat
%    Outputs the two 7x7 matrices corresponding to the 15 wavelets defined  
%    on [0, 15].
%
% M = waveletMat(N)
%    With N an odd integer in [3, 19].
%    Outputs the two (N-1)/2 x (N-1)/2  matrices corresponding to the N  
%    wavelets defined on [0, N].
%
% REFERENCES
%  [1] I. Daubechies and J.C. Lagarias
%        "Two-scale difference equations II. local regularity, infinite
%        products of matrices and fractals",
%      SIAM J. on Math. Anal., 23(4):1031-1079, 1992
%
%  [2] G. Gripenberg, 
%        "Computing the joint spectral radius", 
%      Linear Algebra and its Applications, 234:43-60, 1996
%
%  [3] R. Jungers, 
%        "The Joint Spectral Radius: Theory and Applications" 
%      Vol. 385 section 2.3.3 in Lecture Notes in Control and Information
%      Sciences, Springer-Verlag. Berlin Heidelberg, June 2009
%

load coeffwav

if nargin<1
    N = 15;
else
    if mod(N,2)==0 || N>19 || N<3 || mod(N,1)
        error('The argument should be a nonnegative odd integer in [3, 19]')
    end
end

c = C(1:N+1,(N+1)/2);

T=zeros(N,N,2);
for i=1:N
    for j=1:N
        if (((2*i-j-1)>=0)&&((2*i-j-1)<=N))
            T(i,j,1)=c(2*i-j);
        end
        if (((2*i-j)>=0)&&((2*i-j)<=N))
            T(i,j,2)=c(2*i-j+1);
        end
%         if (((2*i-j-1)>=0)&((2*i-j-1)<=N))
%              Tpsi(i,j,1)=c(N-(2*i-j-1)+1)*(-1)^(N-(2*i-j));
%             %psi est la matrices avec les moins qui permet de determiner
%             %psi quand on connait phi
%         end
%         if (((2*i-j)>=0)&((2*i-j)<=N))
%               Tpsi(i,j,2)=c(N-(2*i-j)+1)*(-1)^(N-(2*i-j+1));
%         end
    end
end

L=(N+1)/2-1;
U=zeros(N,L+1);
for j=1:L+1
    for i=1:N
        U(i,j)=i^(j-1);
    end
end
[q,r]=qr(U);
projecteur=q(:,L+2:N)' ;

Tproj(:,:,1)=projecteur*T(:,:,1)*projecteur';
Tproj(:,:,2)=projecteur*T(:,:,2)*projecteur';

M = [{Tproj(:,:,1)}, {Tproj(:,:,2)}];

end