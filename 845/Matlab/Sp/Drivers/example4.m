loaddata;
clear opt;
opt.initialvec=ones(size(A,1),1);
opt.tolerance=1e-5;
eigenvalues=eigifp(A,B,opt)
quit;
