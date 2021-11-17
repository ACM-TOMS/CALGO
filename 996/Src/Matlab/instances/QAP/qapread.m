function [objPoly, I01, Icomp, penaltyPoly, C, dd] = qapread(instance)
% [objPoly, Icomp, I01, penaltyPoly] = qapread('./qapdata/nug12.dat');
% min trace(XAXB') s.t. X'X = I, X_{ij} \in {0,1}.
%  -> min   x' * kron(B, A) * x    s.t.    x \in {0,1}, .......

  filename = [instance '.dat'];

  fid = fopen(filename,'r');
  if (fid == -1); error('file cannot be opened'); end

  [datavec, ~] = fscanf(fid, '%c');
  fclose('all'); 
  
  linefeeds = strfind(datavec, char(10));
  datavec(linefeeds) = blanks(length(linefeeds)); 
  datavec = sscanf(datavec, '%f'); 

  nn = datavec(1); 
  n2 = nn * nn; 
  A = datavec(2:(n2 + 1)); 
  B = datavec((n2 + 2):(2 * n2 + 1)); 
  
  A = reshape(A, nn, nn);
  B = reshape(B, nn, nn); 
  
  %% objPoly
  Q = kron(B, A);
  objPoly = quad2Poly(Q, [], []);
  
  %% complementarity
  % offdiag(X'X) = 0
  colCompPattern = repmat(speye(nn), nn, nn) - speye(n2);
  rowCompPattern = kron(speye(nn), ones(nn) - speye(nn));
  compPattern = colCompPattern | rowCompPattern;
  [rowIdx, colIdx] = find(triu(compPattern));
  compIdx = repmat(1:length(rowIdx), 2, 1);
  Icomp = sparse(compIdx, [rowIdx'; colIdx'], true, length(rowIdx), n2);
  
  %% binary
  I01 = true(1, n2);
  
  %% diag(X'X) = diag(I).
  C = [kron(ones(1, nn), speye(nn)); kron(speye(nn), ones(1, nn))];
  dd = ones(2 * nn, 1);
  penaltyPoly = quad2Poly(C' * C, -2 * (C' * dd), dd' * dd);
end