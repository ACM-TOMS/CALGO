function [A, b, c, K] = SDPAtoSeDuMi(SDPA)

%
% This function convets SDPA sparse format into SeDuMi Date.
%
% K is the structure with field K.s, K.q, K.l and K.f.
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

b = -SDPA.bVect';
if ~issparse(b)
	b = sparse(b);
end
m = SDPA.mDim;
%%% <--- Set K.f, K.l, K.q and K.s%%%
if isfield(SDPA,'typeCone')
  eq_block_idx = find(SDPA.typeCone == -1);
  lp_block_idx = find(SDPA.typeCone ==  1);
  socp_block_idx = find(SDPA.typeCone ==  2);
  sdp_block_idx = find(SDPA.typeCone == 3);
else
  eq_block_idx = find(SDPA.blockStruct < -1);
  lp_block_idx = find(SDPA.blockStruct == -1);
  sdp_block_idx = find(SDPA.blockStruct > 1);
  socp_block_idx = [];
end
K.f = 0;
K.l = 0;
K.q = [];
K.s = [];

if ~isempty(eq_block_idx) 
  K.f = -sum(SDPA.blockStruct(eq_block_idx),2)/2;
end
if ~isempty(lp_block_idx) 
  K.l =  -sum(SDPA.blockStruct(lp_block_idx),2);
end
if isempty(sdp_block_idx) && isempty(socp_block_idx)
  n = K.f + K.l;
elseif isempty(sdp_block_idx) && ~isempty(socp_block_idx)
  K.q = SDPA.blockStruct(socp_block_idx);
  n = K.f + K.l + K.q*K.q';
elseif ~isempty(sdp_block_idx) && isempty(socp_block_idx)
  K.s = SDPA.blockStruct(sdp_block_idx);
  n = K.f + K.l + K.s*K.s';  
else
  K.s = SDPA.blockStruct(sdp_block_idx);
  K.q = SDPA.blockStruct(socp_block_idx);
  n = K.f + K.l + K.q*K.q' + K.s*K.s';
end
%%% Set K.f, K.l and K.s --->%%%

%%
%% renumber block No. with the order 'equation', 'linear
%% inequality', 'SDP'
%% 
no_eq_block = length(eq_block_idx);
no_lp_block = length(lp_block_idx);
no_socp_block = length(socp_block_idx);
no_sdp_block = length(sdp_block_idx);
matNo(eq_block_idx) = (1:no_eq_block);
matNo(lp_block_idx) = no_eq_block+(1:no_lp_block);
matNo(socp_block_idx) = no_eq_block+no_lp_block+(1: ...
		    no_socp_block);
matNo(sdp_block_idx) = no_eq_block+no_lp_block+no_socp_block+(1: ...
		    no_sdp_block);
        
SDPA.sparseMatrix(:,2) = matNo(SDPA.sparseMatrix(:,2)); 
SDPA.blockStruct = [-SDPA.blockStruct(eq_block_idx)/2,-SDPA.blockStruct(lp_block_idx), SDPA.blockStruct(socp_block_idx),SDPA.blockStruct(sdp_block_idx)];

%%
%% gather the blocks derived from polynomial equation and
%% inequality with basisSupports = (1) into one diagonal block
%%  

SDPAversion = 6.2;
if ~isempty(eq_block_idx) && SDPAversion <= 6.2
  %%
  %% In SDPA 6.2, we need to separate equation into two inequality.
  %% we delete negative part derived by dividing one
  %% equation into two inequalities.
  %%
  EqBlockIdx = find(SDPA.sparseMatrix(:,2) <= no_eq_block);
  EqBlockNo = SDPA.sparseMatrix(EqBlockIdx,2);
  EqBlockSt = SDPA.blockStruct(EqBlockNo)';
  if size(EqBlockSt,1) ~= size(EqBlockIdx,1)
     EqBlockSt = EqBlockSt'; 
  end
  minusPartIdx = find(EqBlockSt < SDPA.sparseMatrix(EqBlockIdx,3));
  SDPA.sparseMatrix(EqBlockIdx(minusPartIdx),:) = [];
end
if  ~isempty(eq_block_idx) || ~isempty(lp_block_idx)
  %%
  %% gather some block derived from polynomial equations or
  %% inequalities into one block.
  %%
  EqLpBlockIdx = find(SDPA.sparseMatrix(:,2) <= no_eq_block+no_lp_block);
  SdpBlockIdx = find(SDPA.sparseMatrix(:,2) > no_eq_block+no_lp_block);
  EqLpBlockNo = SDPA.sparseMatrix(EqLpBlockIdx,2);
  EqLpBlockSt = SDPA.blockStruct(1:no_eq_block+no_lp_block);
  EqLpBlockSt = cumsum(EqLpBlockSt);
  EqLpBlockSt(end) = 0;
  EqLpBlockSt = circshift(EqLpBlockSt',1);
  
  SDPA.sparseMatrix(EqLpBlockIdx,3) =  EqLpBlockSt(EqLpBlockNo) ...
      +SDPA.sparseMatrix(EqLpBlockIdx,3);
  SDPA.sparseMatrix(EqLpBlockIdx,4) = SDPA.sparseMatrix(EqLpBlockIdx,3);
  SDPA.sparseMatrix(EqLpBlockIdx,2) = 1;
  SDPA.sparseMatrix(SdpBlockIdx,2) = 1-no_eq_block-no_lp_block+ ...
      SDPA.sparseMatrix(SdpBlockIdx,2);
  IdxVec = [0;K.q';K.s'];
  
else
  IdxVec = [K.q';K.s'];
end
%%
%% Number all elements in data matrices.
%%
IdxVecSq = cumsum(IdxVec.^2);
IdxVecSq = K.f+K.l+IdxVecSq;
IdxVecSq(end) = 0;
IdxVecSq = circshift(IdxVecSq,1);

%%
%% Make coefficient 'c' of objective function. 
%%
F0_idx = find(SDPA.sparseMatrix(:,1)==0); 
Info_F = SDPA.sparseMatrix(F0_idx,:);
OffdiagIdx = find(Info_F(:,3) ~= Info_F(:,4));
idx = [1,2,4,3,5];
Info_F = [Info_F;Info_F(OffdiagIdx,idx)];

Row = Info_F(:,3);
Col = Info_F(:,4);
Oneidx = (Row-1).*IdxVec(Info_F(:,2)) + Col+IdxVecSq(Info_F(:,2));
c = sparse(Oneidx,1,-Info_F(:,5),n,1);

%%
%% Make coefficient 'A' of constraints.
%%
F_idx = find(SDPA.sparseMatrix(:,1));
Info_F = SDPA.sparseMatrix(F_idx,:);
OffdiagIdx = find(Info_F(:,3) ~= Info_F(:,4));
Info_F = [Info_F;Info_F(OffdiagIdx,idx)];
Row = Info_F(:,3);
Col = Info_F(:,4);
Oneidx = (Row-1).*IdxVec(Info_F(:,2)) + Col+IdxVecSq(Info_F(:,2));
A = sparse(Info_F(:,1),Oneidx,-Info_F(:,5),m,n); 

%%
%% Change some K.s into K.q
%%

if ~isempty(K.q)
   idx = [];
   for i=1:length(K.q)
      idx = [idx, K.q(1:i-1)*K.q(1:i-1)'+(K.q(i)+1:K.q(i)^2)]; 
   end
   idx = K.f+K.l+idx;
   c(idx) = [];
   A(:,idx) = [];
end

return;

% $Header: /home/waki9/CVS_DB/SparsePOPdev/subPrograms/SDPAtoSeDuMi.m,v 1.1.1.1 2007/01/11 11:31:50 waki9 Exp $
