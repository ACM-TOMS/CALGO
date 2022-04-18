function experiments
% we have to make this a function, because subfunctions are not allowed
% in scripts

disp('Part I: One pole, n = 100.');
al = [5, 5i, 1.001, 0.001i]; % four different cases
for j = 1:length(al),
    disp(sprintf('\nPole is %1.3f + %1.3fi',real(al(j)),imag(al(j))));
    a = al(j)*ones(1,100);
    b = 1 ./ (a + (2*(real(a)>=0) - 1) .* sqrt(a.^2 - 1));
    b = b(1:end-1).';
    for w = 1:3,
	disp(sprintf('Weight is %1i',w));
	[x,l,e] = rcheb(a,w);
	disp(sprintf('Maximum error on zeros: %1.3e',max(abs(e))));
	m = max(abs(1 - (f(b,x,w)*l').*(1 - abs(b).^2)/(2*pi)));
	disp(sprintf('Maximum error in quadrature sum: %1.3e',m));
    end
end

disp(sprintf('\nPart II: A few distinct poles, each with multiplicity 100.'));
al = {[2,-2], [i,inf], [1.1, 0.1i, -1.1], [1.1, 0.1i, -1.1, inf]};
for j = 1:length(al),
    l = length(al{j});
    as(1:2:2*l-1) = real(al{j}); as(2:2:2*l) = imag(al{j});
    disp([sprintf('\nPoles are '),sprintf('%1.3f + %1.3fi, ',as)]);
    a = repmat(al{j},1,100);
    b = 1 ./ (a + (2*(real(a)>=0) - 1) .* sqrt(a.^2 - 1));
    b = b(1:end-1).';
    for w = 1:3,
	disp(sprintf('Weight is %1i',w));
	[x,l,e] = rcheb(a,w);
	disp(sprintf('Maximum error on zeros: %1.3e',max(abs(e))));
	m = max(abs(1 - (f(b,x,w)*l').*(1 - abs(b).^2)/(2*pi)));
	disp(sprintf('Maximum error in quadrature sum: %1.3e',m));
    end
end

disp(sprintf('\nPart III: All poles distinct.'));
al = {[1:10]*i*0.001, 1/2*(exp(i*pi*[0:0.01:2])/3 ...
    + 3*exp(-i*pi*[0:0.01:2])), [-1:0.01:1] + i};
as = {'at the integer multiples of 0.001i.', 'on an ellipse.', ...
    'on the horizontal line [-1,1] + i.'};
for j = 1:length(al),
    disp([sprintf('\nPoles are '),as{j}]);
    a = al{j};
    b = 1 ./ (a + (2*(real(a)>=0) - 1) .* sqrt(a.^2 - 1));
    b = b(1:end-1).';
    for w = 1:3,
	disp(sprintf('Weight is %1i',w));
	[x,l,e] = rcheb(a,w);
	disp(sprintf('Maximum error on zeros: %1.3e',max(abs(e))));
	m = max(abs(1 - (f(b,x,w)*l').*(1 - abs(b).^2)/(2*pi)));
	disp(sprintf('Maximum error in quadrature sum: %1.3e',m));
    end
end

disp(sprintf('\nPart IV: Poles extremely close to the boundary.'));
al = {repmat([-0.6:0.2:0.6] + 1e2*eps*i,1,10), 1/2*( ...
    0.9999995*exp(i*pi*[0:0.01:2]) + exp(-i*pi*[0:0.01:2])/0.9999995)};
as = {'on the horizontal line [-0.8,0.8] + 1e2*eps*i', 'on an ellipse very close to [-1,1]'};
for j = 1:length(al),
    disp([sprintf('\nPoles are '),as{j}]);
    a = al{j};
    b = 1 ./ (a + (2*(real(a)>=0) - 1) .* sqrt(a.^2 - 1));
    b = b(1:end-1).';
    for w = 1:3,
	disp(sprintf('Weight is %1i',w));
	[x,l,e] = rcheb(a,w);
	disp(sprintf('Maximum error on zeros: %1.3e',max(abs(e))));
	m = max(abs(1 - (f(b,x,w)*l').*(1 - abs(b).^2)/(2*pi)));
	disp(sprintf('Maximum error in quadrature sum: %1.3e',m));
	disp(sprintf('Number of nodes where full accuracy was not reached: %d (out of %d).', length(find(abs(e) > 50*eps)), length(x)));
    end
end



function fval = f(b,x,w)
% compute the functions abs(phi_k(x))^2, as discussed in the companion
% paper, for k = 1:length(a)-1

bm = repmat(b,1,length(x));
z = exp(i*acos(x)); zm = repmat(z,length(b),1);

B = cumprod((zm(1:end-1,:) - bm(1:end-1,:)) ./ ...
    (1 - conj(bm(1:end-1,:)).*zm(1:end-1,:)));
B = [ones(1,length(z)); B];
Bs = cumprod((zm(1:end-1,:) - conj(bm(1:end-1,:))) ./ ...
    (1 - bm(1:end-1,:).*zm(1:end-1,:)));
Bs = [ones(1,length(z)); Bs];
    
c = (-1)^floor(w/2);

fval = sqrt(2)^(w-1)*(zm.^w.*Bs./(1-bm.*zm) + c./((zm - bm).*B));
if w > 1,
    fval = fval./(zm-1);
    if w > 2,
	fval = fval./(zm+1);
    end
end

fval = fval.*conj(fval);
