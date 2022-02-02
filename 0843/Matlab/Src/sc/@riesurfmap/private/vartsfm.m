function [z,zb] = vartsfm(y,n)

cs = cumsum(cumprod([1;exp(-y(1:n-3))]));
theta = pi*cs(1:n-3)/cs(length(cs));
z = ones(n,1);
z(1:n-3) = exp(i*theta);
z(n-2:n-1) = [-1;-i];
%%ey = exp(y(n-2:end));
%%r = ey(1:2:end);
%%t = ey(2:2:end);
%%%zb = (r+0.05)./(r+1).*exp(i*2*pi*(t)./(t+1));
%%zb = (r-1)./(r+1).*exp(i*pi*t./(t+1));
%%%zb = exp(-r+i*pi*(t-1)./(t+1));

u = y(n-2:2:end);
v = y(n-1:2:end).^2;
zb = (u+i*v-i)./(u+i*v+i);

%%nb = (length(y)-n+3)/2;
%%u = exp(y(n-2:n-2+nb-1)) + i*exp(y(n-2+nb:end));
%%u = u.^2;
%%zb = (u-i)./(u+i);
