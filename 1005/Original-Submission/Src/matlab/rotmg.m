function [d, x1, param] = rotmg(d, x)
  % Translation of a version of srotmg.f available on calgo.acm.org. This
  % earlier version is free of the bug that the netlib versions have (LAPACK v.
  % 3.2 through 3.8.0).

  [d1, d2] = deal(d(1), d(2));
  [x1, y1] = deal(x(1), x(2));
  gam = 4096;
  gamsq = gam^2;
  rgamsq = 1/gamsq;
  
  if d1 < 0, zero_h_d_and_sx1(); return, end
  % case-d1-nonnegative
  p2 = d2*y1;
  if p2 == 0
    flag = -2;
    param(1) = flag;
    return
  end
  % regular-case..
  p1 = d1*x1;
  q2 = p2*y1;
  q1 = p1*x1;
  if abs(q1) > abs(q2)
    h12 = p2/p1;
    h21 = -y1/x1;
    u = 1 - h12*h21;
    if u <= 0, zero_h_d_and_sx1(); return, end
    flag = 0;
    d1 = d1/u;
    d2 = d2/u;
    x1 = x1*u;
  else
    if q2 < 0, zero_h_d_and_sx1(); return, end
    flag = 1;
    h11 = p1/p2;
    h22 = x1/y1;
    u = 1 + h11*h22;
    temp = d2/u;
    d2 = d1/u;
    d1 = temp;
    x1 = y1*u;
  end
  scale_check()
  finish()
  return

  function finish
    if flag < 0
      param(2) = h11;
      param(3) = h21;
      param(4) = h12;
      param(5) = h22;
    elseif flag == 0
      param(3) = h21;
      param(4) = h12;
    else
      param(2) = h11;
      param(5) = h22;
    end
    param(1) = flag;
    d = [d1; d2];
  end
  
  function scale_check
    %gam1 = 1;
    %gam2 = 1;
    %flg = flag;
    while d1 <= rgamsq
      if d1 == 0, break, end
      fix_h()
      d1 = d1*gam^2;
      x1 = x1/gam;
      h11 = h11/gam;
      h12 = h12/gam;
      %gam1 = gam1*gam;
    end
    while d1 >= gamsq
      fix_h()
      d1 = d1/gam^2;
      x1 = x1*gam;
      h11 = h11*gam;
      h12 = h12*gam;
      %gam1 = gam1*gam;
    end
    while abs(d2) <= rgamsq
      if abs(d2) == 0, break, end
      fix_h()
      d2 = d2*gam^2;
      h21 = h21/gam;
      h22 = h22/gam;
      %gam2 = gam2*gam;
    end
    while abs(d2) >= gamsq
      fix_h()
      d2 = d2/gam^2;
      h21 = h21*gam;
      h22 = h22*gam;
      %gam2 = gam2*gam;
    end
    %fprintf('flag, gam1, gam2 = %.0f, %.0f, %.0f\n', flg, gam1, gam2)
    finish()
  end
  
  function zero_h_d_and_sx1
    flag = -1;
    h11 = 0;
    h12 = 0;
    h21 = 0;
    h22 = 0;
    d1 = 0;
    d2 = 0;
    x1 = 0;
    finish()
  end
  
  function fix_h
    if flag == 0
      h11 = 1;
      h22 = 1;
    elseif flag > 0
      h21 = -1;
      h12 = 1;
    end
    flag = -1;
  end
end
