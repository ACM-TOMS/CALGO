% TEST_DER2ARRAY_MDS_SET  Check functions der2array and mds_set

function test_der2array_mds_set
  fprintf('TESTING DER2ARRAY AND MDS_SET.... ');
  test(0,1);
  test(1,1);
  test(0,2);
  test(1,2);
  test(2,2);
  disp('OK')
end

function test(ppq,r);
  nPar = (ppq)*r^2 + r*(r+1)/2;
  F = rand(r,r);  
  Fd = rand(r,r,nPar);
  % PARTIAL TEST OF MDS_SET
  FS = mds_set(F,Fd);  
  ascertain(almostequal(FS.mat,F));
  ascertain(almostequal(FS.der{1}{1,1},Fd(:,:,1)));
  if ppq>0, ascertain(almostequal(FS.der{2}{1,1},Fd(:,:,r^2+1))); end
  if ppq>0 && r>1
    ascertain(almostequal(FS.der{1}{2,1},reshape(Fd(:,:,2),r,r)));
    ascertain(almostequal(FS.der{1}{1,2},reshape(Fd(:,:,r+1),r,r)));
  end
  % TEST DER2ARRAY FOR ONE-ELEMENT INPUT
  [F1,Fd1] = der2array(FS);
  Fd2 = der2array(FS);
  ascertain(almostequal(F,F1{1}));
  ascertain(almostequal(Fd,Fd1{1}) && almostequal(Fd,Fd2{1}));
  % TEST DER2ARRAY FOR TWO-ELEMENT INPUT
  G = rand(r,r);
  Gd = rand(r,r,nPar);
  GS = mds_set(G,Gd);
  [F1,Fd1] = der2array([FS GS]);
  ascertain(almostequal(F,F1{1}));
  ascertain(almostequal(G,F1{2}));
  ascertain(almostequal(Fd,Fd1{1})); 
  ascertain(almostequal(Gd,Fd1{2}));
end
