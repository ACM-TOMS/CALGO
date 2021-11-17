package body Complex_Laurent_Polynomial_Systems is 

-- CREATORS :

  function Create ( p : Laur_Sys ) return Eval_Laur_Sys is
    
    res : Eval_Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Create(p(k));
    end loop;
    return res;
  end Create;

  function Create ( p : Laur_Sys ) return Eval_Coeff_Laur_Sys is
    
    res : Eval_Coeff_Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Create(p(k));
    end loop;
    return res;
  end Create;

-- ADDITIONAL ARITHMETIC OPERATIONS :

  function "*" ( a : double_complex; p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := a*p(k);
    end loop;
    return res;
  end "*";

  function "*" ( p : Laur_Sys; a : double_complex ) return Laur_Sys is
  begin
    return a*p;
  end "*";

  procedure Mult_Cmplx ( p : in out Laur_Sys; a : in double_complex ) is
  begin
    for k in p'range loop
      Mult_Coeff(p(k),a);
    end loop;
  end Mult_Cmplx;

  function Eval ( p : Laur_Sys; x : double_complex; i : natural )
                return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for j in p'range loop
      res(j) := Eval(p(j),x,i);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Laur_Sys; x : Vector ) return Vector is

    res : Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Eval_Laur_Sys; x : Vector ) return Vector is

    res : Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Eval(p(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( p : Eval_Coeff_Laur_Sys;
                  c : Complex_Vectors_of_Vectors.Vector; x : Vector )
                return Vector is

    res : Vector(p'range);

  begin
    for i in p'range loop
      res(i) := Eval(p(i),c(i).all,x);
    end loop;
    return res;
  end Eval;

  function Diff ( p : Laur_Sys; i : natural ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for j in p'range loop
      res(j) := Diff(p(j),i);
    end loop;
    return res;
  end Diff;

  procedure Diff ( p : in out Laur_Sys; i : in natural ) is
  begin
    for j in p'range loop
      Diff(p(j),i);
    end loop;
  end Diff;

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Laur_Sys ) is
  begin
    for k in p'range loop
      Clear(p(k));
    end loop;
  end Clear;

  procedure Clear ( p : in out Eval_Coeff_Laur_Sys ) is
  begin
    for k in p'range loop
      Clear(p(k));
    end loop;
  end Clear;

end Complex_Laurent_Polynomial_Systems;
