function y=Horner(coef,n,t)

coef_aux=coef;

y=coef_aux(n+1)*ones(size(t));
for r=n:-1:1
    y=y.*t+coef_aux(r);
end
