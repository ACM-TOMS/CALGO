with Random_Number_Generators;           use Random_Number_Generators;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers;
with Float_Vectors_of_Vectors;

package body Float_Lifting_Functions is

-- AUXILIARIES :

  function Flt2Cmplx ( x : FLoat_Vectors.Vector )
                     return Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with complex entries.

    res : Complex_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := CMPLX(x(i));
    end loop;
    return res;
  end Flt2Cmplx;

-- RANDOM FLOATING-POINT LIFTING :

  function Random_Lift ( lflow,lfupp : double_float ) return double_float is

    res : double_float := random;                          -- in [-1,1]

  begin
    res := ((1.0+res)/2.0)*lflow + ((1.0-res)/2.0)*lfupp;  -- in [lflow,lfupp]
    return res;
  end Random_Lift;

  function Random_Lift ( v : Vector; lflow,lfupp : double_float )
                       return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v; 
    res(res'last) := Random_Lift(lflow,lfupp);
    return res;
  end Random_Lift;

  function Random_Lift ( l : List; lflow,lfupp : double_float ) return List is

    res,res_last : List;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Random_Lift(Head_Of(tmp).all,lflow,lfupp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Random_Lift;

  function Random_Lift ( l : Arrays_of_Float_Vector_Lists.Array_of_Lists;
                         lflow,lfupp : Vector )
                       return Arrays_of_Float_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Float_Vector_Lists.Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Random_Lift(l(i),lflow(i),lfupp(i));
    end loop;
    return res;
  end Random_Lift;

-- LINEAR LIFTING FUNCTIONS :

  function Linear_Lift ( x,v : Vector ) return Vector is

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := x*v;
    return res;
  end Linear_Lift;

  function Linear_Lift ( f : Face; v : Vector ) return Face is

    res : Face := new Float_Vectors_of_Vectors.Vector(f'range);

  begin
    for i in res'range loop
      res(i) := new Float_Vectors.Vector'(Linear_Lift(f(i).all,v));
    end loop;
    return res;
  end Linear_Lift;

  function Linear_Lift ( l : List; v : Vector ) return List is

  -- DESCRIPTION :
  --   Returns a linearly lifted list of points.

    res,res_last : List;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Linear_Lift(Head_Of(tmp).all,v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Linear_Lift;

  function Linear_Lift ( f : Faces; v : Vector ) return Faces is

    res,res_last : Faces;
    tmp : Faces := f;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Linear_Lift(Head_Of(tmp),v));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Linear_Lift;

-- RANDOM FLOATING-POINT LINEAR LIFTING :

  function Random ( n : natural; lflow,lfupp : double_float ) return Vector is

    res : Vector(1..n);

  begin
    for i in res'range loop
      res(i) := Random_Lift(lflow,lfupp);
    end loop;
    return res;
  end Random;

-- POLYNOMIAL LIFTING FUNCTIONS :

  function Polynomial_Lift ( lf : Poly; x : Vector ) return Vector is

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := REAL_PART(Eval(lf,Flt2Cmplx(x)));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Poly; x : Vector ) return Vector is

    res : Vector(x'first..x'last+1);

  begin
    res(x'range) := x;
    res(res'last) := REAL_PART(Eval(lf,Flt2Cmplx(x)));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Poly; l : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Poly; l : List ) return List is

    res,res_last,tmp : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),l(i));
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Polynomial_Lift(lf(i),l(i));
    end loop;
    return res;
  end Polynomial_Lift;

end Float_Lifting_Functions;
