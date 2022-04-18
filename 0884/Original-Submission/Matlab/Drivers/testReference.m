% Test for reference.m
%
% Given an arbitrary triangle, the Argyris basis is computed
% analitically using the symbolic toolbox. 
%
% Next we evalute numerically and analitically the basis at 
% several random points and compare the results.

p=rand(2,3);  % an arbitrary triangle

[m,argy]=reference(p);
argyx =diff(argy,'x');
argyy =diff(argy,'y');
argyxx=diff(argyx,'x');
argyxy=diff(argyx,'y');
argyyy=diff(argyy,'y');
argy  =vectorize(inline(argy));
argyx =vectorize(inline(argyx));
argyy =vectorize(inline(argyy));
argyxx=vectorize(inline(argyxx));
argyxy=vectorize(inline(argyxy));
argyyy=vectorize(inline(argyyy));

% 20 arbitrary points on the triangle
q=rand(20,3); 
q=inv(diag(sum(q,2)))*q;  % barycentric coordinates
z1=p(1,:)*q'; 
z2=p(2,:)*q';
  
aux=evalArgy(z1,z2,p);
aux2=argy(z1,z2);
disp('maximum of error')
disp(max(abs(aux(:)-aux2(:))))

[auxx,auxy]=evalGradArgy(z1,z2,p);

aux2x=argyx(z1,z2); aux2y=argyy(z1,z2);
disp('maximum of the error for the first derivatives')
disp(max(max(abs([auxx(:)-aux2x(:) auxy(:)-aux2y(:)]))))

[auxxx auxxy auxyy]=evalHessArgy(z1,z2,p);
aux2xx=argyxx(z1,z2); aux2xy=argyxy(z1,z2); aux2yy=argyyy(z1,z2);

disp('maximum of the error for the second derivatives')
disp(max(max(abs([auxxx(:)-aux2xx(:) auxxy(:)-aux2xy(:) auxyy(:)-aux2yy(:)]))))
