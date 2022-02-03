with Greatest_Common_Divisor;           use Greatest_Common_Divisor;
with Float_Linear_Inequality_Solvers;   use Float_Linear_Inequality_Solvers;

--with text_io,integer_io;                use text_io,integer_io;
--with Integer_Vectors_io;                use Integer_Vectors_io;
--with Float_Vectors_io;                  use Float_Vectors_io;

package body Float_to_Integer_Rounding is

-- BASIC CONVERSION OPERATIONS :

  function Round ( v : Float_Vectors.Vector ) return Integer_Vectors.Vector is

    res : Integer_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := integer(v(i));
    end loop;
    return res;
  end Round;

  function Convert_to_Float ( v : Integer_Vectors.Vector )
                            return Float_Vectors.Vector is

    res : Float_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := double_float(v(i));
    end loop;
    return res;
  end Convert_to_Float;

  function Convert_to_Float ( m : Integer_Matrices.matrix )
                            return Float_Matrices.matrix is

  -- DESCRIPTION :
  --   Converts an integer matrix to a floating point matrix.

    res : Float_Matrices.matrix(m'range(1),m'range(2));

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        res(i,j) := double_float(m(i,j));
      end loop;
    end loop;
    return res;
  end Convert_to_Float;

-- AUXILIARIES FOR CONTINUED FRACTIONS :

  function Truncate ( x : double_float ) return integer is

  -- DESCRIPTION :
  --   Truncates x to the nearest lower integer number.

    res : integer := integer(x);

  begin
    if double_float(res) <= x
     then return res;
     else return res-1;
    end if;
  end Truncate;

  function lcm ( v : Integer_Vectors.Vector ) return integer is

  -- DESCRIPTION :
  --   Returns the least common multiple of all entries in the vector.

  -- REQUIRED : for all i in v'range: v(i) /= 0.

    res : integer := v(v'first);

  begin
    for i in v'first+1..v'last loop
      res := lcm(res,v(i));
    end loop;
    return res;
  end lcm;

  procedure Swap ( i1,i2 : in out integer ) is

  -- DESCRIPTION : Swaps the integers i1 and i2.

    tmp : integer := i1;

  begin
    i1 := i2; i2 := tmp;
  end Swap;

  procedure Positive_Fractions
                 ( x,tol : in double_float; a,b : out integer ) is

  -- DESCRIPTION :
  --   By continued fractions, a and b are computed such that
  --   abs(x - a/b) < tol.

  -- REQUIRED : x > 0.

    xx : double_float := x;
    al : integer := Truncate(xx);

  begin
    if abs(xx - double_float(al)) < tol
     then a := al; b := 1;
     else xx := 1.0/(xx - double_float(al));
          declare
            g1,g2,h1,h2 : integer;
          begin
            g1 := 1; g2 := 0; h1 := 0; h2 := 1;
            g2 := al*g1 + g2;
            h2 := al*h1 + h2;
            Swap(g1,g2); Swap(h1,h2);
            loop
              al := Truncate(xx);
              g2 := al*g1 + g2;
              h2 := al*h1 + h2;
              Swap(g1,g2); Swap(h1,h2);
              exit when (abs(xx - double_float(al)) < tol);
              exit when (abs(x - double_float(g1)/double_float(h1)) < tol);
              xx := 1.0/(xx - double_float(al));
            end loop;
            a := g1; b := h1;
          end;
    end if;
  end Positive_Fractions;

-- CONVERSIONS BY CONTINUED FRACTIONS :

  procedure Fractions ( x,tol : in double_float; a,b : out integer ) is
  begin
    if x >= 0.0
     then Positive_Fractions(x,tol,a,b);
     else declare
            aa : integer;
          begin
            Positive_Fractions(-x,tol,aa,b);
            a := -aa;
          end;
    end if;
  end Fractions;

  function Scale_to_Integer ( v : Float_Vectors.Vector; tol : in double_float )
                            return Integer_Vectors.Vector is

    res,den : Integer_Vectors.Vector(v'range);
    comden : integer;

  begin
   -- put("The floating point solution : "); put(v,3,3,3); new_line;
    for i in v'range loop
      Fractions(v(i),tol,res(i),den(i));
    end loop;
   -- put("The nominators : "); put(res); new_line;
   -- put("The denominators : "); put(den); new_line;
    comden := lcm(den);
   -- put("The common denominator : "); put(comden,1); new_line;
    if comden /= 1
     then for i in res'range loop
            res(i) := res(i)*(comden/den(i));
          end loop;
    end if;
    return res;
  exception
    when numeric_error => return (res'range => 0);
  end Scale_to_Integer;

  function Scale_to_Integer
               ( ine : Float_Matrices.Matrix; v : Float_Vectors.Vector;
                 tol : in double_float ) return Integer_Vectors.Vector is

    res : Integer_Vectors.Vector(v'range);
    fltres : Float_Vectors.Vector(res'range);
    inctol : double_float := 0.1;

  begin
    res := Round(v);
    fltres := Convert_to_Float(res);
    while not Satisfies(ine,fltres,tol) loop
      res := Scale_to_Integer(v,inctol);
      exit when (inctol <= tol);
      fltres := Convert_to_Float(res);
      inctol := inctol*0.1;
    end loop;
    return res;
  end Scale_to_Integer;

end Float_to_Integer_Rounding;
