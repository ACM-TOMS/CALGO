function showStruct(sadata, varargin)
%showStruct displays the structure of the DAE in the current figure, or a 
%new figure if none is open.
%showStruct(sadata)
%showStruct(sadata, options)
%sadata is an object returned by daeSA.
%
%This function can have one to three optional key-value pair arguments,
%where the keys are 'disptype', 'submat', and 'blocksubmat'.
%
%[1] 'disptype','original'   displays the original DAE structure (default)
%    'disptype','blocks'     displays coarse block-triangularized (BTS)
%    'disptype','fineblocks' displays fine BTS
%
%[2]  'submat', submat
%Corresponding to 'submat' is a vector submat of size 2 or 4 specifying
%what submatrix to be displayed.
%If submat is of length 2, diagonal submatrix comprising rows submat(1)
%to submat(2) and columns submat(1) to submat(2) is displayed.
%If submat is of length 4, submatrix comprising rows submat(1) to
%submat(2) and columns submat(3) to submat(4) is displayed.
%
%Example:
%   submat = [1 3];
%   showStruct(sadata,'submat',submat);
%   showStruct(sadata,'disptype','blocks','submat',submat);
%   showStruct(sadata,'disptype','fineblocks','submat',submat);
%
%[3] 'blocksubmat',blocksubmat
%Similarly, corresponding to 'blocksubmat' is a vector of size 2 or 4
%specifying which blocks to be displayed.
%If blocksubmat is of length 2, diagonal blocks comprising block rows
%blocksubmat(1) to blocksubmat(2) and block columns blocksubmat(1) to
%blocksubmat(2) is displayed.
%If blocksubmat is of length 4, diagonal blocks comprising block rows
%blocksubmat(1) to blocksubmat(2) and block columns submat(3) to submat(4)
%is displayed.
%
%This pair argument MUST be used with displaying coarse/fine BTS.
%Example:
%   blocksubmat=[1 3]    
%   showStruct(sadata,'disptype','blocks','blocksubmat',blocksubmat);
%   showStruct(sadata,'disptype','fineblocks','blocksubmat',blocksubmat);
%
%The pairs 'submat',submat and 'blocksubmat',blocksubmat CANNOT
%both exist in the argument list 'options'.
%
%If the DAE is not structurally well posed (SWP) due to missing
%equations and/or variables, they are displayed with reds dots in the
%figure.  If the DAE has more equations than the argument n passed to
%daeSA, then the DAE is not SWP, and no figure is dislpayed.  If it
%has more variables than n, then Matlab reports an error.
%
%See also daeSA, isSWP, getOffsets, getBTF
%
%Copyright 2012 Guangning Tan, Ned Nedialkov, and John Pryce
DAESAshowMAIN(sadata, varargin{:})
end