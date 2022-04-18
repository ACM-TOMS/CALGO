function krontests 
%% Table 1
disp('Table 1. Timing comparisons with full Kronecker product matrix')
disp(' d,	m=n,  reps,       tf,       tb,       tr,       tk,  (2tk)/(tf+tb)')
repvec=[30 20 10];
vmax=0;
for d=[2 3];
 for i=[10 20 30]
   m=repmat(i,1,d);
   n=m;
   reps=repvec(i/10);
   [A,v]=gendata(m,n);
   [tt,vv]=krontest1(A,v,reps);
   if max(vv)>vmax, vmax=max(vv); end
   if any(vv>1e-10)
     fprintf('results differ by more than 10^-10')
   end
   tt(1:4)=tt(1:4)/reps;
   % displays timing results
   fprintf('%2.0f,  %2.0f,    %2.0f, %8.4f, %8.4f, %8.4f, %8.4f,    %8.4f \n', [d i reps tt])  
 end
end
fprintf('maximum relative difference in output: %1.4e\n',vmax)

%% Table 2
disp('Table 2. Timing comparisons with kronm & kronmult - square equally sized matrices')
disp(' d,	m=n,  reps,       tf,       tb,       tm,       tt,  (2tm)/(tf+tb)')
repvec=[20 15 10; 12 8 4; 6 4 2];
vmax=0;
for d=[3 4 5];
 for i=[30 40 50]
   m=repmat(i,1,d);
   n=m;
   reps=repvec(d-2,i/10-2);
   [A,v]=gendata(m,n);
   [tt,vv]=krontest2(A,v,reps);
   if max(vv)>vmax, vmax=max(vv); end
   if any(vv>1e-10)
     fprintf('results differ by more than 10^-10')
   end
   tt(1:4)=tt(1:4)/reps;
   % displays timing results
   fprintf('%2.0f,  %2.0f,    %2.0f, %8.4f, %8.4f, %8.4f, %8.4f,       %8.4f \n', [d i reps tt])  
 end
end
fprintf('maximum  relative difference in output: %1.4e\n',vmax)

%% Table 3
disp('Table 3. Timing comparisons with kronm - rectangular unequally sized matrices')
disp('forward approach is optimal')
disp('d,  reps,	    tf,       tb,       tm,  max(tf,tb)/min(tf,tb), tm/min(tf,tb)')
mm=[10 20 30 40 50];
nn=[50 40 30 20 10];
repvec=[30 20 10];
vmax=0;
for d=[3 4 5];
 m=mm(end-d+1:end);
 n=nn(1:d);
 reps=repvec(d-2);
 [A,v]=gendata(m,n);
 [tt,vv]=krontest3(A,v,reps);
 if max(vv)>vmax, vmax=max(vv); end
 if any(vv>1e-10)
   fprintf('results differ by more than 10^-10')
 end
 tt(1:3)=tt(1:3)/reps;
 % displays timing results
 fprintf('%2.0f,   %2.0f, %8.4f, %8.4f, %8.4f,              %8.4f,       %8.4f \n', [d reps tt])     
end
 
disp('backward approach is optimal')
disp('d,  reps,	    tf,       tb,       tm, max(tf,tb)/min(tf,tb), tm/min(tf,tb)')
for d=[3 4 5];
 n=mm(end-d+1:end);
 m=nn(1:d);
 reps=repvec(d-2);
 [A,v]=gendata(m,n);
 [tt,vv]=krontest3(A,v,reps);
 if max(vv)>vmax, vmax=max(vv); end
 if any(vv>1e-10)
   fprintf('results differ by more than 10^-10')
 end
 tt(1:3)=tt(1:3)/reps;
 % displays timing results
fprintf('%2.0f,   %2.0f, %8.4f, %8.4f, %8.4f,              %8.4f,       %8.4f \n', [d reps tt])     
end
fprintf('maximum  relative difference in output: %1.4e\n',vmax) 

%%
 
% generates a set of matrices with matrix i m(i) x n(i)
 function [A,v]=gendata(m,n)
   d=length(m); 
   A=cell(1,d); 
   for i=1:d; 
     A{i}=randn(m(i),n(i)); 
   end
   v=randn(prod(n),1);

 % performs the tests for table 1
 function [tt,vv]=krontest1(A,v,reps)
 fopt=struct('forward',true);  
 bopt=struct('forward',false);  
 t1=0; t2=0; t3=0; t4=0;
 for i=1:reps
   tic; v1=ckronx(A,v,fopt);   t1=t1+toc;   % my algorithm - forward
   tic; v2=ckronx(A,v,bopt);   t2=t2+toc;   % my algorithm - backward
   tic; AA=ckron(A);           t3=t3+toc;   % form A
   tic; v3=AA*v;               t4=t4+toc;   % A*v
 end
 tt=[t1 t2 t3 t4 t4/mean([t1 t2])];
 % test to ensure results are the same (except for round off)
 errcheck=@(v1,v2) max(abs(v1(:)-v2(:))./max(abs(v1(:))));
 vv=[errcheck(v1,v2) errcheck(v1,v3) errcheck(v2,v3)];
 

% performs the tests for table 2 
function [tt,vv]=krontest2(A,v,reps)
 fopt=struct('forward',true);  
 bopt=struct('forward',false);  
 Ar=A(end:-1:1);                % kronm performs (Ad x ... x A1)*B
  
 t1=0; t2=0; t3=0; t4=0;
 for i=1:reps
   tic; v1=ckronx(A,v,fopt);  t1=t1+toc;   % my algorithm - forward
   tic; v2=ckronx(A,v,bopt);  t2=t2+toc;   % my algorithm - backward
   tic; v3=kronm(Ar,v);       t3=t3+toc;   % kronm algorithm from matlab central
   tic; v4=kronmult(A,v);     t4=t4+toc;   % kronmult algorithm from matlab central
 end
 tt=[t1 t2 t3 t4 t3/mean([t1 t2])];
 % test to ensure results are the same (except for round off)
 errcheck=@(v1,v2) max(abs(v1(:)-v2(:))./max(abs(v1(:))));
 vv=[errcheck(v1,v2) errcheck(v1,v3) errcheck(v1,v4)];


 % performs the tests for table 3 
 function [tt,vv]=krontest3(A,v,reps)
 fopt=struct('forward',true);  
 bopt=struct('forward',false);  
 Ar=A(end:-1:1);                % kronm performs (Ad x ... x A1)*B
  
 t1=0; t2=0; t3=0;
 for i=1:reps
   tic; v1=ckronx(A,v,fopt);  t1=t1+toc;   % my algorithm - forward
   tic; v2=ckronx(A,v,bopt); t2=t2+toc;   % my algorithm - backward
   tic; v3=kronm(Ar,v);         t3=t3+toc;   % kronm algorithm from matlab central
 end
 tt=[t1 t2 t3 max(t1,t2)/min(t1,t2) t3/min(t1,t2)];
 % test to ensure results are the same (except for round off)
 errcheck=@(v1,v2) max(abs(v1(:)-v2(:))./max(abs(v1(:))));
 vv=[errcheck(v1,v2) errcheck(v1,v3) errcheck(v2,v3)];
 
 % computes a sequential Kronecker product 
 function AA=ckron(A)
   AA=A{1};
   for i=2:length(A);
     AA=kron(AA,A{i});
   end

