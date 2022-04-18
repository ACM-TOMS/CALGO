function [da, xa, info] = rotmg_rmd(d, x1, param, da, x1a, parama)
  % Reverse mode derivative of rotmg. 
  %
  % Input parameters:
  %    d, x1, param:  The results returned by the corresponding rotmg call
  %    da:            Adjoint of the d returned by rotmg
  %    x1a:           Adjoint of the x1 returned by rotmg
  %    parama:        Adjoint of param in entries that contain H-elements
  %                   according to the value of flag = param(1)
  %
  % Output parameters:
  %    da:            Adjoint of the d supplied to rotmg (assigned to)
  %    xa             Adjoint of the x supplied to rotmg (assigned to)
  %    info           Information on detected parameters:
  %                     info.flag: Flag before scaling
  %                     info.gam1: Gamma for d1
  %                     info.gam2: Gamma for d2
  %                   (Note: d1 and d2 are swapped when Flag = 1)
  
  flag = param(1);
  [d1, d2] = deal(d(1), d(2));
  [d1a, d2a] = deal(da(1), da(2));
  [gam1, gam2] = deal(1);
  [d1u, d2u] = deal(d1, d2);
  if flag == -2
    xa = [x1a; 0];
  elseif flag == -1
    [h11, h21, h12, h22] = deal(param(2), param(3), param(4), param(5));
    [h11a,h21a,h12a,h22a] = deal(parama(2),parama(3),parama(4),parama(5));
    % Use D = H'*D~*H and H*[x1;y1] = [x1~;0] (where D~ and x1~ are values on
    % exit from rotmg to reconstruct D, x1 and y1 as they were on entry to rotmg
    % Also reconstruct the scales, gam1 and gam2, applied to d1 and d2 by rotmg
    % if such scaling was done.
    u = h11*h22 - h12*h21; % det(H)
    x1r = h22*x1/u;
    y1r = -h21*x1/u;
    d1r = h11^2*d1 + h21^2*d2;
    d2r = h12^2*d1 + h22^2*d2;
    % Unscale:
    if  abs(d1r*x1r^2) > abs(d2r*y1r^2)
      flag = 0;
      gam1 = 1/h11;
      gam2 = 1/h22;
      h12 = h12*gam1;
      h21 = h21*gam2;
      h12a = h12a/gam1;
      h21a = h21a/gam2;
      u = 1 - h12*h21;
    else
      flag = 1;
      gam1 = 1/h12;
      gam2 = -1/h21;
      h11 = h11*gam1;
      h22 = h22*gam2;
      h11a = h11a/gam1;
      h22a = h22a/gam2;
      u = 1 + h11*h22;
    end
    %x1 = x1*gam1;
    x1a = x1a/gam1;
    d1u = d1/gam1^2;
    d2u = d2/gam2^2;
    d1a = d1a*gam1^2;
    d2a = d2a*gam2^2;
  elseif flag == 0
    [h21, h12] = deal(param(3), param(4));
    [h21a, h12a] = deal(parama(3), parama(4));
    u = 1 - h12*h21;
    % Reconstruct d1, d2, x1 and y1 as they were on entry to rotmg:
    d1r = d1*u;
    d2r = d2*u;
    x1r = x1/u;
    y1r = -h21*x1r;
  elseif flag == 1
    [h11, h22] = deal(param(2), param(5));
    [h11a, h22a] = deal(parama(2), parama(5));
    u = 1 + h11*h22;
    % Reconstruct d1, d2, x1 and y1 as they were on entry to rotmg:
    d1r = d2*u;
    d2r = d1*u;
    y1r = x1/u;
    x1r = h22*y1r;    
  end
  info = struct('gam1', gam1, 'gam2', gam2, 'flag', flag);
  %print()
  if flag == 0
    p1 = d1r*x1r;
    ua = x1r*x1a;
    x1a = u*x1a;
    d2a = d2a/u;
    ua = ua - d2u*d2a;
    d1a = d1a/u;
    ua = ua - d1u*d1a;
    h12a = h12a - h21*ua;
    h21a = h21a - h12*ua;
    y1a = -h21a/x1r;
    x1a = x1a + h21*y1a;
    p2a = h12a/p1;
    p1a = -h12*p2a;
  elseif flag == 1
    p2 = d2r*y1r;
    y1a = u*x1a;
    ua = y1r*x1a;
    temp = d1u;
    tempa = d1a;
    d1a = d2a/u;
    ua = ua - d2u*d1a;
    d2a = tempa/u;
    ua = ua - temp*d2a;
    h11a = h11a + h22*ua;
    h22a = h22a + h11*ua;
    x1a = h22a/y1r;
    y1a = y1a - h22*x1a;
    p1a = h11a/p2;
    p2a = -h11*p1a;
  elseif flag == -2
    return
  else
    error('Illegal value of flag');
  end
  d1a = d1a + x1r*p1a;
  x1a = x1a + d1r*p1a;
  d2a = d2a + y1r*p2a;
  y1a = y1a + d2r*p2a;
  da = [d1a d2a];
  xa = [x1a y1a];

  function print() %#ok<DEFNU>
    fprintf('\nIn rotmg_rmd\nReconstructed variables at rotmg entry:\n');
    fprintf('  d = %.6f, %.6f\n', d1r, d2r);
    fprintf('  x1r = %.6f, y1r = %.6f\n', x1r, y1r);
    fprintf('Reconstructed variables before rotmg scaling:\n');
    fprintf('  flag = %d\n', flag);
    fprintf('  d = %.6f, %.6f\n', d1u, d2u);
    fprintf('  x1 = %.6f\n', x1);
    if flag == 0
      fprintf('  H = [1 %.4f; %.4f 1]\n', h12, h21);
    else
      fprintf('  H = [%.4f 1; -1 %.4f]\n', h11, h22);
    end
    if flag == 0
      fprintf('  Ha = [1 %.6f; %.6f 1]\n', h12a, h21a);
    else
      fprintf('  Ha = [%.6f 1; -1 %.6f]\n', h11a, h22a);
    end
    fprintf('Adjoints after unscaling:\n');
    fprintf('  da = %.6f, %.6f\n', d1a, d2a);
    fprintf('  x1a = %.6f\n', x1a);
  end
end
