%%********************************************************
%%
%% [Q] = biqread('~/SDPdata/BIQ/be100.1.sparse'); 
%%
%%********************************************************

  function [Q] = biqread(instance)

  fid = fopen(instance,'r');
  if (fid == -1); error('file cannot be opened'); end

  [datavec,count] = fscanf(fid,'%c');
  fclose(fid); 
  linefeeds = findstr(datavec,char(10));
  datavec(linefeeds) = blanks(length(linefeeds)); 
  datavec = sscanf(datavec,'%f'); 

  n = datavec(1); nz = datavec(2); 
  tmp = datavec(3:end);   
  ii = tmp(1:3:end-2);
  jj = tmp(2:3:end-1); 
  vv = tmp(3:3:end); 
  Q = spconvert([ii,jj,vv;n,n,0]); 
  Q = triu(Q,1) + triu(Q,1)' + spdiags(diag(Q),0,n,n); 
%%********************************************************
