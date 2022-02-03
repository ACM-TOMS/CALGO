with Floating_Point_Numbers;    use Floating_Point_Numbers;
with Complex_Numbers;           use Complex_Numbers;

package Random_Number_Generators is

  function Random ( lower, upper : integer ) return integer;

  -- DESCRIPTION :
  --   Returns an integer between the lower and upper bound,
  --   bounds also included, randomly generated.

  function Random return double_float;

  -- DESCRIPTION :
  --   Returns a randomly chosen floating point number
  --   between -1.0 and 1.0.

  function Random return double_complex;

  -- DESCRIPTION :
  --   Using plainly the floating point random number generator
  --   for generating real and imaginary part,
  --   a randomly chosen complex number is constructed.

  function Random ( modulus : double_float ) return double_complex;

  -- DESCRIPTION :
  --   Generates a random complex number with a given modulus,
  --   so only the argument angle will be chosen at random.

  function Random1 return double_complex;

  -- DESCRIPTION :
  --   Generates a random complex number with modulus one.

end Random_Number_Generators;
