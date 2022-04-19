loaddata;
eigenvalues=eigifp(A)
eigenvalues=eigifp(A,B)
eigenvalues=eigifp(A,B,3)
clear opt;
opt.initialvec=ones(size(A,1),1);
opt.tolerance=1e-5;
eigenvalues=eigifp(A,B,opt)
opt.useprecon='NO';
eigenvalues=eigifp(A,B,opt)
opt.useprecon=0;
eigenvalues=eigifp(A,B,opt)
clear opt;
opt.inneriteration=32;
eigenvalues=eigifp(A,B,opt)
quit;
