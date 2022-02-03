with Mathematical_Functions;         use Mathematical_Functions;
with Integer_Vectors;
with Complex_Linear_System_Solvers;  use Complex_Linear_System_Solvers;  

package body Predictors is

-- PREDICTORS for t :

  procedure Real_Predictor
              ( t : in out double_complex; target : in double_complex;
                h,tol : in double_float; pow : in positive := 1;
                hh : out double_float ) is

    nt : double_complex;

  begin
    if pow = 1
     then nt := t + CMPLX(h);
     else nt := CMPLX((h+REAL_PART(t)**pow)**(1.0/double_float(pow)));
    end if;
    if REAL_PART(nt) >= REAL_PART(target)
      or else abs( REAL_PART(nt) - REAL_PART(target) ) <= tol
     then hh := REAL_PART(target) - REAL_PART(t);
          t := target;
     else hh := h;
          t := nt;
    end if;
  end Real_Predictor;

  procedure Complex_Predictor
               ( t : in out double_complex; target : in double_complex;
                 h,tol : in double_float; hh : out double_float;
                 distance : in double_float; trial : in natural ) is

    nt,target_direction,step : double_complex;
    dre,dim,d,x,y,alfa : double_float;
    tr3 : natural range 0..2;

  begin
    target_direction := target - t;
    d := modulus(target_direction);
    if d < tol
     then nt := target - CMPLX(distance);
          hh := modulus(nt - t);
     else tr3 := trial mod 3;
	  case tr3 is
	    when 0 => target_direction := target_direction / CMPLX(d);
                      if (d - h) > distance
                       then step := CMPLX(h) * target_direction;
                       else step := (d - distance) * target_direction;
                      end if;
                      nt := t + step;
	    when 1 => dre := REAL_PART(target_direction);
                      dim := IMAG_PART(target_direction);
                      alfa := ARCTAN(dim/dre);
                      x := h * COS(alfa + PI/4.0);
                      y := h * SIN(alfa + PI/4.0);
                      step := CMPLX(x,y);
                      nt := t + step;
                      if (modulus(nt) > modulus(target)) 
                         or else (modulus(nt - target) < distance)
                        then step := CMPLX(0.0,y) * CMPLX(h);
                             nt := t + step;
                      end if;
	    when 2 => dre := REAL_PART(target_direction);
                      dim := IMAG_PART(target_direction);
                      alfa := ARCTAN(dim/dre);
                      x := h * COS(alfa - PI/4.0);
                      y := h * SIN(alfa - PI/4.0);
                      step := CMPLX(x,y);
                      nt := t + step;
                      if (modulus(nt) > modulus(target)) 
                         or else (modulus(nt - target) < distance)
                        then step := CMPLX(0.0,-y) * CMPLX(h);
                             nt := t + step;
                      end if;
          end case;
          hh := modulus(step);
    end if;
    t := nt;
  end Complex_Predictor;

  procedure Circular_Predictor
              ( t : in out double_complex; theta : in out double_float;
                t0_min_target,target : in double_complex; 
                h : in double_float ) is

    e_i_theta : double_complex;

  begin
    theta := theta + h;
    e_i_theta := CMPLX(COS(theta),SIN(theta));
    t := target + t0_min_target * e_i_theta;
  end Circular_Predictor;

  procedure Geometric_Predictor
              ( t : in out double_complex; target : in double_complex;
                h,tol : in double_float ) is

    nt : double_complex;

  begin
    nt := target - CMPLX(h)*(target - t);
    if abs( REAL_PART(nt) - REAL_PART(target) ) <= tol
     then t := target;
     else t := nt;
    end if;
  end Geometric_Predictor;

-- PREDICTORS FOR x :

  procedure Secant_Predictor
              ( x : in out Solution_Array; prev_x : in Solution_Array;
                fac : in double_complex; dist_x : in double_float ) is

    j : natural;
    xx : Solution_Array(x'range);

  begin
    Copy(x,xx);
    for i in x'range loop
      x(i).v := x(i).v + fac * ( x(i).v - prev_x(i).v );
      j := xx'first;
      Equals(xx,x(i).v,i,dist_x,j);
      if j /= i
       then Copy(xx,x);
            exit;
      end if;
    end loop;
    Clear(xx);
  end Secant_Predictor;

-- PREDICTORS FOR t AND x :

  procedure Secant_Single_Real_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out double_complex; prev_t,target : in double_complex;
                 h,tol : in double_float; pow : in positive := 1 ) is

    hh,tmp : double_float;
    factor : double_complex;

  begin
    tmp := modulus(t - prev_t);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol
     then factor := CMPLX(hh/tmp);
          x := x + factor * ( x - prev_x );
    end if;
  end Secant_Single_Real_Predictor;

  procedure Secant_Multiple_Real_Predictor
               ( x : in out Solution_Array; prev_x : in Solution_Array;
                 t : in out double_complex; prev_t,target : in double_complex;
                 h,tol,dist_x : in double_float; pow : in positive := 1 ) is

    hh,tmp : double_float;
    factor : double_complex;

  begin
    tmp := modulus(t - prev_t);
    Real_Predictor(t,target,h,tol,pow,hh);
    if tmp > tol
     then factor := CMPLX(hh/tmp);
          Secant_Predictor(x,prev_x,factor,dist_x);
    end if;
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Secant_Multiple_Real_Predictor;

  procedure Secant_Single_Complex_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out double_complex; prev_t,target : in double_complex;
                 h,tol,dist_t : in double_float; trial : in natural ) is

    hh,tmp : double_float;
    factor : double_complex;

  begin
    tmp := modulus(t - prev_t);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    if tmp > tol
     then factor := CMPLX(hh/tmp);
          x := x + factor * ( x - prev_x );
    end if;
  end Secant_Single_Complex_Predictor;

  procedure Secant_Multiple_Complex_Predictor
               ( x : in out Solution_Array; prev_x : in Solution_Array;
                 t : in out double_complex; prev_t,target : in double_complex;
                 h,tol,dist_x,dist_t : in double_float; trial : in natural ) is

    hh,tmp : double_float;
    factor : double_complex;

  begin
    tmp := modulus(t - prev_t);
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    if tmp > tol
     then factor := CMPLX(hh/tmp);
          Secant_Predictor(x,prev_x,factor,dist_x);
    end if;
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Secant_Multiple_Complex_Predictor;

  procedure Secant_Circular_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out double_complex; theta : in out double_float;
                 prev_t,t0_min_target,target : in double_complex;
                 h,tol : in double_float ) is

    nt,factor : double_complex;
    tmp : double_float := modulus(t-prev_t);

  begin
    if tmp <= tol
     then Circular_Predictor(t,theta,t0_min_target,target,h);
     else factor := CMPLX(h/tmp);
          Circular_Predictor(t,theta,t0_min_target,target,h);
          x := x + factor * ( x - prev_x );
    end if;
  end Secant_Circular_Predictor;

  procedure Secant_Geometric_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out double_complex; prev_t,target : in double_complex;
                 h,tol : in double_float ) is

    dist_prev,dist : double_float;
    factor,tmp : double_complex;

  begin
    dist_prev := modulus(t - prev_t);     -- old stepsize
    tmp := t;
    Geometric_Predictor(t,target,h,tol);
    dist := modulus(t - tmp);             -- new stepsize
    if dist_prev > tol
     then factor := CMPLX(dist/dist_prev);
          x := x + factor * ( x - prev_x );
    end if;
  end Secant_Geometric_Predictor;

  procedure Tangent_Single_Real_Predictor
               ( x : in out Vector; t : in out double_complex;
                 target : in double_complex; h,tol : in double_float;
                 pow : in positive := 1 ) is

    n : natural := x'last;
    info : integer;
    norm_tan,hh : double_float;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    prev_t : double_complex := t;

  begin
    Real_Predictor(t,target,h,tol,pow,hh);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0
     then rhs := dH(x,prev_t);
          Min_Vector(rhs);
          lusolve(j,n,ipvt,rhs);
          norm_tan := Norm(rhs);
          if norm_tan > tol
           then hh := hh / norm_tan;
                Mult_Coeff(rhs,CMPLX(hh));
                Plus_Vector(x,rhs);
          end if;
    end if;
  end Tangent_Single_Real_Predictor;

  procedure Tangent_Multiple_Real_Predictor
               ( x : in out Solution_Array; t : in out double_complex;
                 target : in double_complex; h,tol,dist_x : in double_float;
                 nsys : in out natural; pow : in positive := 1 ) is

    n : natural := x(x'first).n;
    norm_tan,hh : double_float;
    rhs : Vector(1..n);
    info : integer;
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    xx : Solution_Array(x'range);
    jj : natural;
    prev_t : double_complex := t;

  begin
    Real_Predictor(t,target,h,tol,pow,hh);
    Copy(x,xx);
    for i in x'range loop
      j := dH(x(i).v,prev_t);
      lufac(j,n,ipvt,info);
      if info = 0
       then rhs := dH(x(i).v,prev_t);
            Min_Vector(rhs);
            lusolve(j,n,ipvt,rhs);
            nsys := nsys + 1;
            norm_tan := Norm(rhs);
            if norm_tan > tol
             then hh := hh / norm_tan;
                  Mult_Coeff(rhs,CMPLX(hh));
                  Plus_Vector(x(i).v,rhs);
            end if;
            jj := xx'first;
            Equals(xx,x(i).v,i,dist_x,jj);
            if jj /= i
             then Copy(xx,x); exit;
            end if;
       else Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Tangent_Multiple_Real_Predictor;

  procedure Tangent_Single_Complex_Predictor 
               ( x : in out Vector; t : in out double_complex;
                 target : in double_complex;
                 h,tol,dist_t : in double_float; trial : in natural) is

    n : natural := x'last;
    info : integer;
    norm_tan,hh : double_float;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    prev_t : double_complex := t;

  begin
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0
     then rhs := dH(x,prev_t);
          Min_Vector(rhs);
          lusolve(j,n,ipvt,rhs);
          norm_tan := Norm(rhs);
          if norm_tan > tol
           then hh := hh / norm_tan;
                Mult_Coeff(rhs,CMPLX(hh));
                Plus_Vector(x,rhs);
          end if;
    end if;
  end Tangent_Single_Complex_Predictor;

  procedure Tangent_Multiple_Complex_Predictor
               ( x : in out Solution_Array; t : in out double_complex;
                 target : in double_complex;
                 h,tol,dist_x,dist_t : in double_float;
                 trial : in natural; nsys : in out natural) is

    n : natural := x(x'first).n;
    norm_tan,hh : double_float;
    rhs : Vector(1..n);
    info : integer;
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    xx : Solution_Array(x'range);
    jj : natural;
    prev_t : double_complex := t;

  begin
    Complex_Predictor(t,target,h,tol,hh,dist_t,trial);
    Copy(x,xx);
    for i in x'range loop
      j := dH(x(i).v,prev_t);
      lufac(j,n,ipvt,info);
      if info = 0
       then rhs := dH(x(i).v,prev_t);
            Min_Vector(rhs);
            lusolve(j,n,ipvt,rhs);
            nsys := nsys + 1;
            norm_tan := Norm(rhs);
            if norm_tan > tol
             then hh := hh / norm_tan;
                  Mult_Coeff(rhs,CMPLX(hh));
                  Plus_Vector(x(i).v,rhs);
            end if;
            jj := xx'first;
            Equals(xx,x(i).v,i,dist_x,jj);
            if jj /= i
             then Copy(xx,x); exit;
            end if;
       else Copy(xx,x); exit;
      end if;
    end loop;
    Clear(xx);
    for k in x'range loop
      x(k).t := t;
    end loop;
  end Tangent_Multiple_Complex_Predictor;

  procedure Tangent_Circular_Predictor 
              ( x : in out Vector; t : in out double_complex;
                theta : in out double_float;
                t0_min_target,target : in double_complex;
                h,tol : in double_float ) is

    n : natural := x'last;
    info : integer;
    norm_tan : double_float;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);

  begin
    lufac(j,n,ipvt,info);
    if info = 0
     then rhs := dH(x,t);
          Min_Vector(rhs);
          lusolve(j,n,ipvt,rhs);
    end if;
    Circular_Predictor(t,theta,t0_min_target,target,h);
    if info = 0
     then norm_tan := Norm(rhs);
          if norm_tan > tol
           then Mult_Coeff(rhs,CMPLX(h/norm_tan));
                Plus_Vector(x,rhs);
          end if;
    end if;
  end Tangent_Circular_Predictor;

  procedure Tangent_Geometric_Predictor
               ( x : in out Vector; t : in out double_complex;
                 target : in double_complex; h,tol : in double_float ) is

    n : natural := x'last;
    info : integer;
    norm_tan,step : double_float;
    rhs : Vector(x'range);
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    prev_t : double_complex := t;

  begin
    Geometric_Predictor(t,target,h,tol);
    j := dH(x,prev_t);
    lufac(j,n,ipvt,info);
    if info = 0
     then rhs := dH(x,prev_t);
          Min_Vector(rhs);
          lusolve(j,n,ipvt,rhs);
          norm_tan := Norm(rhs);
          if norm_tan > tol
           then step := modulus(t-prev_t) / norm_tan;
                Mult_Coeff(rhs,CMPLX(step));
                Plus_Vector(x,rhs);
          end if;
    end if;
  end Tangent_Geometric_Predictor;

  function Hermite ( t0,t1,t,x0,x1,v0,v1 : double_complex )
                   return double_complex is

  -- DESCRIPTION :
  --   Returns the value of the third degree interpolating polynomial x(t),
  --   such that x(t0) = x0, x(t1) = x1, x'(t0) = v0 and x'(t1) = v1.

  -- REQUIRED : t0 /= t1.

  -- IMPLEMENTATION :
  --   x(t) = a3*t^3 + a2*t^2 + a1*t + a0,
  --   The four interpolation conditions lead to a linear system in
  --   the coefficients of x(t).  This system is first solved explicitly
  --   and then the polynomial x(t) is evaluated. 

    a0,a1,a2,a3,t10,v10 : double_complex;

  begin
    t10 := t1 - t0;
    v10 := (x1 - x0)/t10;
    a3 := (v1 + v0 - CMPLX(2.0)*v10)/(t10**2);
    a2 := (v10 - v0 - (t1**2 + t1*t0 - CMPLX(2.0)*t0**2)*a3)/t10;
    a1 := v0 - (CMPLX(3.0)*a3*t0 + CMPLX(2.0)*a2)*t0;
    a0 := x0 - ((a3*t0 + a2)*t0 + a1)*t0;
    return (((a3*t + a2)*t + a1)*t + a0);
  end Hermite;

  function Hermite ( t0,t1,t : double_complex; x0,x1,v0,v1 : Vector )
                   return Vector is

  -- DESCRIPTION :
  --   Returns the value of the third degree interpolating polynomial x(t),
  --   such that x(t0) = x0, x(t1) = x1, x'(t0) = v0 and x'(t1) = v1,
  --   for every component.

  -- REQUIRED : t0 /= t1.

    res : Vector(x0'range);

  begin
    for i in res'range loop
      res(i) := Hermite(t0,t1,t,x0(i),x1(i),v0(i),v1(i));
    end loop;
    return res;
  end Hermite;

  procedure Hermite_Single_Real_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out double_complex; prev_t,target : in double_complex;
                  v : in out Vector; prev_v : in Vector;
                  h,tol : in double_float; pow : in positive := 1 ) is

    n : natural := x'last;
    info : integer;
    hh : double_float;
    j : Matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    t1 : double_complex := t;

  begin
    Real_Predictor(t,target,h,tol,pow,hh);
    j := dH(x,t1);
    lufac(j,n,ipvt,info);
    if info = 0
     then v := dH(x,t1);
          Min_Vector(v);
          lusolve(j,n,ipvt,v);
          if abs(prev_t - t1) > tol
           then x := Hermite(prev_t,t1,t,prev_x,x,prev_v,v);
          end if;
    end if;
  end Hermite_Single_Real_Predictor;

end Predictors;
