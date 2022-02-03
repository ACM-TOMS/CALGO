function [z,bkerr,cond] = rtsinfo(p)

% roots, backward error and condition number using "roots"

n = length(p)-1;
z = [roots(p),ones(n,1)];
p = p/p(1);
cond = spcond(z);

[yy,jj]=sort(abs(z)); y = z(jj);
ff = [1]; w = ones(1,n+1);
for k = 1:n
    ff = conv(ff,poly(y(k)));
    if abs(p(k+1)) > 1 
        w(k+1) = 1/abs(p(k+1));
    end
end

bkerr = norm( (ff-p).*w, inf);
