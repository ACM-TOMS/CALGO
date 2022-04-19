with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers;
with Integer_Vectors_Utilities;          use Integer_Vectors_Utilities;
with Random_Number_Generators;

package body Integer_Lifting_Functions is

-- AUXILIARIES :

  function Convert ( v : Integer_Vectors.Vector )
                   return Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the given vector into a vector with complex entries.

    res : Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := CMPLX(double_float(v(i)));
    end loop;
    return res;
  end Convert;

  function Random_Vector ( m : natural; low,upp : integer ) return Vector is

  -- DESCRIPTION :
  --   Returns a random vector of range 1..m of randomly generated integer
  --   values between low and upp.

    res : Vector(1..m);

  begin
    for k in res'range loop
      res(k) := Random_Number_Generators.Random(low,upp);
    end loop;
    return res;
  end Random_Vector;

  function Random ( vec : Integer_Vectors.Vector ) return integer is

  -- DESCRIPTION :
  --   Returns a random number from the given vector.

    index : integer := Random_Number_Generators.Random(vec'first,vec'last);

  begin
    return vec(index);
  end Random;

-- LINEAR LIFTING :

  function Linear_Lift ( lf,v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := lf*v;
    return res;
  end Linear_Lift;

  function Linear_Lift ( lf : Vector; l : List ) return List is

    res,res_last : List;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Linear_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Linear_Lift;

  function Linear_Lift ( lf : Integer_Vectors_of_Vectors.Vector;
                         l : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Linear_Lift(lf(i).all,l(i));
    end loop;
    return res;
  end Linear_Lift;   

-- POLYNOMIAL LIFTING :

  function Polynomial_Lift ( lf : Poly; v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := integer(REAL_PART(Eval(lf,Convert(v))));
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Poly; v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := integer(REAL_PART(Eval(lf,Convert(v))));
    return res;
  end Polynomial_Lift;
 
  function Polynomial_Lift ( lf : Poly; l : List ) return List is
 
    res : List;
    elf : Eval_Poly := Create(lf);

  begin
    res := Polynomial_Lift(elf,l);
    Clear(elf);
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Eval_Poly; l : List ) return List is

    res,res_last : List;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Polynomial_Lift(lf,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Polynomial_Lift;

  function Polynomial_Lift ( lf : Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists is

    res : Array_of_Lists(l'range);
    elf : Eval_Poly_Sys(lf'range) := Create(lf);
  
  begin
    res := Polynomial_Lift(elf,l);
    Clear(elf);
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

-- RANDOM LIFTING :

  function Random_Lift ( lflow,lfupp : integer; v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := Random_Number_Generators.Random(lflow,lfupp);
    return res;
  end Random_Lift;

  function Random_Lift ( lflow,lfupp : integer; l : List ) return List is

    res,res_last : List;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Random_Lift(lflow,lfupp,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Random_Lift;

  function Random_Lift ( lflow,lfupp : Vector; l : Array_of_Lists )
                       return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Random_Lift(lflow(i),lfupp(i),l(i));
    end loop;
    return res;
  end Random_Lift;

-- RANDOM LINEAR LIFTING :

  function Random_Linear_Lift ( lflow,lfupp : integer; v : Vector )
                              return Vector is

    lf : Vector(v'range) := Random_Vector(v'last,lflow,lfupp);
 
  begin
    return Linear_Lift(lf,v);
  end Random_Linear_Lift;

  function Random_Linear_Lift ( lflow,lfupp : integer; l : List )
                              return List is
  begin
    if Is_Null(l)
     then return l;
     else declare
            n : constant natural := Head_Of(l)'length;
            lf : Vector(Head_Of(l)'range) := Random_Vector(n,lflow,lfupp);
          begin
            return Linear_Lift(lf,l);
          end;
    end if;
  end Random_Linear_Lift;

  function Random_Linear_Lift ( lflow,lfupp : Vector; l : Array_of_Lists )
                              return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Random_Linear_Lift(lflow(i),lfupp(i),l(i));
    end loop;
    return res;
  end Random_Linear_Lift;

-- POINT-WISE LIFTING :

  function Point_Lift ( lf : integer; v : Vector ) return Vector is

    res : Vector(v'first..v'last+1);

  begin
    res(v'range) := v;
    res(res'last) := lf;
    return res;
  end Point_Lift;

  function Point_Lift ( lf : Vector; l : List ) return List is

    res,res_last : List;
    tmp : List := l;
    ind : integer := lf'first;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Point_Lift(lf(ind),Head_Of(tmp).all));
      ind := ind + 1;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Point_Lift;

  function Point_Lift ( lf : Integer_Vectors_of_Vectors.Vector;
                        l : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(l'range);

  begin
    for i in res'range loop
      res(i) := Point_Lift(lf(i).all,l(i));
    end loop;
    return res;
  end Point_Lift;

end Integer_Lifting_Functions;
