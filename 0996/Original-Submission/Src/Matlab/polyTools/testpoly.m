poly1.supports = [1 0 0; 1 1 0; 0 0 1];
poly1.coef = [1 2 3]';
poly2.supports = [1 1 0; 0 0 1; 0 1 0];
poly2.coef = [1 2 3]';

newPoly = addPoly(poly1,poly2);

%%

poly.supports = [1 0 0 0 0; 2 0 0 0 0; 1 1 0 1 1 ];
poly.coef = [1; 1; 1];

newPoly = simplifyPoly(poly,true(1,5),1);

%%
poly.supports = [0 0; 1 0];
poly.coef = [-1; 1];
mulpoly = poly;
for ii=1:1
    mulpoly = multiplyPoly(mulpoly,poly);
end
%%
poly.supports = [1 0; 0 1];
poly.coef = [-1; 1];
newPoly = powPoly(poly,10);

%%
poly.supports = [1 1 1; 1 2 1];
poly.coef = [1;2];
newPoly = varTransPoly(poly,2,2,-1);