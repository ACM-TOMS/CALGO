function y = asech(x);

% AD implementation of asech.m
% Code written by Ilyssa Sanders and Anil V. Rao
% January 2009

y.value = asech(x.value);
outerDerivative = -1./(full(x.value(:)).*(1-full(x.value).^2));
y = compositeDerivative(x,y,outerDerivative);
y = class(y,'ad');
