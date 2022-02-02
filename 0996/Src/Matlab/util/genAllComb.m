function params = genAllComb(paramset)
%GENALLCOMB Generate all combinations of parameters 
%
% Usage: 
%   params = GENALLCOMB(paramset);
% 
% Example:
%   If 
%     paramset{1} = {0, 1};
%     paramset{2} = {'a', 'b'};
%   
%   then
%     params = {0, 'a';...
%             0, 'b';...
%             1, 'a';...
%             1, 'b'};


totalcomb = 1;
for ii = 1:length(paramset)
  totalcomb = totalcomb * length(paramset{ii});
end

radix = zeros(1, length(paramset));
radix(end) = 1;
for ii = length(paramset)-1:-1:1
    radix(ii) = radix(ii+1) * length(paramset{ii+1});
end
totalcomb = radix(ii) * length(paramset{ii});

params = cell(totalcomb, length(paramset));
for ii = 0:totalcomb-1
    tmp = ii;
    for jj = 1:length(paramset)
        idx = floor(ii/radix(jj));
        params{tmp+1,jj} = paramset{jj}{idx+1};
        ii = rem(ii,radix(jj));
    end
end

end

