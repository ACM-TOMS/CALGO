%% Shepp-Logan Phantom in 2D
% For reference see <incontrastsRef.html>.

function q = shepp2D(X1,X2,varargin)
% Shepp-Logan Phantom

%P = phantom('Modified Shepp-Logan',200);
%imshow(P)

% seti = varargin{end};
% R = seti.rCD/2;

n = sqrt(length(X1));

% -- image (square)

%p = phantom('Modified Shepp-Logan',n); % p is n x n (this is not sparse enough)
p = phantom('Modified Shepp-Logan',floor(n/2)); % p is n/2 x n/2

% build mask
x1 = 1:n;
indx = (floor(n/2)-n/4<x1)&(x1<=floor(n/2)+n/4);
mask = zeros(n,n);
mask(indx,indx) = 1;

% use mask
y = zeros(size(mask));
y(mask~=0) = p; % phantom

% -- image end

q = reshape(y,[1 n^2]);

q = double(q);

end

