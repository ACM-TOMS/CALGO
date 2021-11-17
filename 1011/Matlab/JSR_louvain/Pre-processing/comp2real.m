function Mreal = comp2real(M)
%
% MREAL = COMP2REAL(M)
%
% Transforms a set of square nxn matrices into
% a set of (2*n)x(2*n) real matrices.
% 
% MREAL{i} = [A -B;B A];
%
% with A = real(M{i})
%      B = imag(M{i})
  
m = length(M);
Mreal = cell(1,m);

for imat = 1:m
    A = real(M{imat});
    B = imag(M{imat});
    
    Mreal{imat} = [A -B;B A];
        
end


end