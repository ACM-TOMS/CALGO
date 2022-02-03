with Mathematical_Functions;  use Mathematical_Functions;
with Machines;                use Machines;

package body Random_Number_Generators is

  a : constant integer := 13849;
  m : constant integer := 65536;
  c : constant integer := 56963;
  seed : integer := Process_Id;

  function Random ( lower, upper : integer ) return integer is

    f : double_float;

  begin
    f := Random; -- f in [-1,1]
    f := 0.5*(double_float(upper-lower)*f + double_float(lower+upper));
                                         -- f in [lower,upper]
    return integer(f); -- rounding the result to integer number
  end Random;

  function Random return double_float is

    x : double_float;

  begin
    seed := (a*seed + c) mod m;
    x := double_float(seed)/double_float(m);
    x := 2.0 * x - 1.0;
    return x;
  end Random;
  
  function Random return double_complex is
  begin
    return CMPLX(Random,Random);
  end Random;

  function Random ( modulus : double_float ) return double_complex is

    arg : double_float;

  begin
    arg := PI*Random;  -- in [-pi,+pi]
    return CMPLX(modulus*COS(arg),modulus*SIN(arg));
  end Random;

  function Random1 return double_complex is

    arg : double_float; 

  begin
    arg := PI*Random;  -- in [-pi,+pi]
    return CMPLX(COS(arg),SIN(arg));
  end Random1;

end Random_Number_Generators;
