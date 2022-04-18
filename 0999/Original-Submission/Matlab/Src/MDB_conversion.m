function ccd = MDB_conversion(MPd, Hd, MPs, Hs, ccs)

% Conversion from source to destination MDB-spline form;
% the number of segments must be the same in both forms
%
% INPUT
%   MPd   : destination MDB-spline multi-patch
%   Hd    : destination extraction matrix
%   MPs   : source MDB-spline multi-patch
%   Hs    : source extraction matrix
%   ccs   : source coefficient vector
%
% OUTPUT
%   ccd   : destination coefficient vector

m = length(MPd.P);
dds = reshape(ccs, 1, []) * Hs;
ddd = zeros(1, MPd.mu(end));
for i = 1:m
   dds_loc = dds(MPs.mu(i)+1:MPs.mu(i+1));
   ddd_loc = B_conversion(MPd.P(i), MPs.P(i), dds_loc);
   ddd(MPd.mu(i)+1:MPd.mu(i+1)) = ddd_loc;
end
ccd = ddd / Hd;

end
