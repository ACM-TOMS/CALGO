function M = overlapMat
% 
% M = overlapMat   generates set of two 20x20 matrices whose joint spectral 
%                  radius characterizes the asymptotic behaviour of the 
%                  number of so-called overlap-free words of length n, when
%                  n tends to infinity [1].
%
% The spectral radius of this set has been proven to be 
%    sqrt(rho(M{1}*M{2})) = 2.51793404...   [2]
%
% REFERENCES
% [1] R.M. Jungers, V. Protasov, V.D. Blondel,
%      "Overlap-free words and spectra of matrices",
%     Theoretical Computer Science, 410(38):3670-3684, 2009
%
% [2] N. Guglielmi, V. Protasov,
%      "Exact computation of joint spectral characteristics of linear
%      operators",
%     Foundations of Computational Mathematics, 13(1):37-97, 2013

AA=zeros(15,15);%AA,BB,CC sont les matrices A,B,C de cassaigne
BB=zeros(15,15);
CC=zeros(15,15);
AA(1,8)=1;
AA(1,9)=2;
AA(1,10)=1;
AA(2,3)=1;
AA(2,4)=1;
AA(2,6)=1;
AA(2,7)=1;
AA(3,8)=1;
AA(3,9)=1;
AA(5,1)=1;
AA(5,2)=2;
AA(5,5)=1;
AA(6,3)=1;
AA(6,6)=1;
AA(8,8)=1;

BB(1,8)=1;
BB(1,9)=2;
BB(1,10)=1;
BB(3,6)=1;
BB(3,7)=1;
BB(3,14)=1;
BB(3,15)=1;
BB(4,3)=1;
BB(4,4)=1;
BB(4,12)=1;
BB(4,13)=1;
BB(8,5)=1;
BB(9,2)=1;
BB(9,11)=1;
BB(10,1)=1;

CC(1,8)=2;
CC(1,9)=4;
CC(1,10)=2;
CC(2,3)=1;
CC(2,4)=1;
CC(2,6)=1;
CC(2,7)=1;
CC(3,6)=1;
CC(3,7)=1;
CC(3,8)=1;
CC(3,9)=1;
CC(4,3)=1;
CC(4,4)=1;
CC(6,2)=1;
CC(6,5)=1;
CC(7,1)=1;
CC(7,2)=1;
CC(8,6)=2;
CC(8,14)=2;
CC(9,3)=1;
CC(9,12)=1;
CC(11,15)=1;
CC(11,12)=1;
CC(11,13)=1;
CC(11,14)=1;
CC(12,15)=1;
CC(12,14)=1;
CC(13,13)=1;
CC(13,12)=1;
CC(14,11)=1;
CC(15,11)=1;

DD=zeros(15,15);

M{1}=[CC,DD;AA,BB];
M{2}=[AA,BB;DD,CC];

A = zeros(30,30,2);
A(:,:,1)=[CC,DD;AA,BB];
A(:,:,2)=[AA,BB;DD,CC];

M{1}=[A(1:10,1:10,1) A(1:10,16:25,1);A(16:25,1:10,1) A(16:25,16:25,1)];

M{2}=[A(1:10,1:10,2) A(1:10,16:25,2);A(16:25,1:10,2) A(16:25,16:25,2)];

end