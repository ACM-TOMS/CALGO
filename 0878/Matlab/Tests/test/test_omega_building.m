%TEST_OMEGA_BUILDING Test the functions omega_build and omega_remove_miss
%
%  TEST_OMEGA_BUILDING prints out OK (three times) if omega_build works ok for
%  three test cases, all with r=2, p=3, n=7 but q is set to 2, 3 and 4
%  respectively. It may also be useful for debugging.
%
%  TEST_OMEGA_BUILDING QUIET runs more quietly
%
%  TEST_OMEGA_BUILDING MISS tests also omega_remove_miss

function test_omega_building(varargin)
  [varargin,quiet] = getflags(varargin,'quiet');
  [varargin,miss] = getflags(varargin,'miss');
  S0=[3 4; 8 8];
  S1=[3 4; 1 1];
  S2=[3 4; 2 2];
  S3=[3 4; 3 3];
  G0=[4 5; 8 8];
  G1=[4 5; 1 1];
  G2=[4 5; 2 2];
  G3=[4 5; 3 3];
  G4=[4 5; 4 4];
  W0=[5 6; 8 8];
  W1=[5 6; 1 1];
  W2=[5 6; 2 2];
  W3=[5 6; 3 3];
  W4=[5 6; 4 4];
  S={S0 S1 S2 S3};
  G={G0 G1 G2};
  W={W0 W1 W2};
  r=2; p=3; q=2; n = 7;
  if ~miss, fprintf('TESTING OMEGA_BUILD...');
  else      fprintf('TESTING OMEGA_BUILD AND OMEGA_REMOVE_MISS...'); end    
  fprintf_if(~quiet && ~miss, '\n  test q<p...');
  [Su,Olow] = omega_build(S, G, W, p, n);
  O = zeros(r);
  SuOK = [
    S0 O  O
    S1 S0 O
    S2 S1 S0];
  OlowOK = [
    G2 G1 W0
    G2 W1 W0
    W2 W1 W0
    W2 W1 W0];
  ascertain(isequal(tril(SuOK), tril(Su)))
  ascertain(isequal(OlowOK, Olow))
  
  fprintf_if(~quiet && ~miss, 'OK\n  test q=p+1...');
  SOKB = [
    S0 O  O   O
    S1 S0 O   O
    S2 S1 S0  O
    G3 G2 G1 W0];
  GB={G0 G1 G2 G3 G4};
  WB={W0 W1 W2 W3 W4};
  [SuB,OlowB] = omega_build(S, GB, WB, p, n);
  OlowOKB = [
    G4 G3 G2 W1 W0
    G4 G3 W2 W1 W0
    G4 W3 W2 W1 W0];
  ascertain(isequal(tril(SOKB), tril(SuB)))
  ascertain(isequal(OlowOKB, OlowB))

  fprintf_if(~quiet && ~miss, 'OK\n  test q=p+2...');
  SC = {S0 S1 S2};
  [SuC,OlowC] = omega_build(SC, GB, WB, 2, n);
  SOKC = [
    S0  O  O  O
    S1 S0  O  O
    G2 G1 W0  O
    G3 G2 W1 W0];
  OlowOKC = [
    G4 G3 W2 W1 W0
    G4 W3 W2 W1 W0
    W4 W3 W2 W1 W0];
  ascertain(isequal(tril(SOKC), tril(SuC)))
  ascertain(isequal(OlowOK, Olow))

  if miss
    miss = false(r, n);
    miss([1 3 4 7 9 10 14]) = true;
    SuoOK = [
      8 0 0
      4 3 4
      2 8 8];
    OlowoOK = [
      1 1 8
      6 5 6
      2 8 8
      5 6 5];
    ko = find_missing_info(miss);
    [Suo, Olowo] = omega_remove_miss(Su, Olow, miss);
    ascertain(isequal(tril(SuoOK), tril(Suo)));
    ascertain(isequal(OlowoOK, Olowo(:,1:3)))
  end
  disp('OK')
end
