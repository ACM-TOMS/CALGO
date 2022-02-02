with Floating_Point_Numbers;      use Floating_Point_Numbers;

package Complex_Numbers is

  type double_complex is private;
  I: constant double_complex;

-- real => double_complex

  function CMPLX ( X : double_float ) return double_complex; -- along real axis

-- real * real => double_complex :

  function CMPLX ( X,Y : double_float ) return double_complex; -- Cartesian

-- real * argument => double_complex :

  function CMPLX_POLAR ( R : double_float; A : double_float )
                       return double_complex;                  -- Polar

-- double_complex => real :

  function "ABS"   ( X : double_complex ) return double_float;
  function MODULUS ( X : double_complex )
                   return double_float renames "ABS";
  function RADIUS  ( X : double_complex )
                   return double_float renames "ABS";

  function REAL_PART ( X : double_complex ) return double_float;
  function IMAG_PART ( X : double_complex ) return double_float;

-- double_complex => argument :

  function ARGUMENT     ( X : double_complex ) return double_float;
  function ANGLE	( X : double_complex ) 
                        return double_float renames ARGUMENT;

--  double_complex * real => double_complex :

  function "+"  ( X : double_complex; Y : double_float ) return double_complex;
  function "-"  ( X : double_complex; Y : double_float ) return double_complex;
  function "*"  ( X : double_complex; Y : double_float ) return double_complex;
  function "/"  ( X : double_complex; Y : double_float ) return double_complex;

-- inline only for double_complex "/" real => double_complex
-- other "/" overloadings not inline

  pragma inline ("/");

-- real * double_complex => double_complex :

  function "+"  ( X : double_float; Y : double_complex ) return double_complex;
  function "-"  ( X : double_float; Y : double_complex ) return double_complex;
  function "*"  ( X : double_float; Y : double_complex ) return double_complex;
  function "/"  ( X : double_float; Y : double_complex ) return double_complex;

-- double_complex => double_complex :

  function "+"   ( X : double_complex ) return double_complex;
  function "-"   ( X : double_complex ) return double_complex;

-- double_complex * double_complex => double_complex :

  function "+"  ( X,Y : double_complex ) return double_complex;
  function "-"  ( X,Y : double_complex ) return double_complex;
  function "*"  ( X,Y : double_complex ) return double_complex;
  function "/"  ( X,Y : double_complex ) return double_complex;

-- inline for all overloadings of named functions :

  pragma inline(CMPLX, "ABS", MODULUS, ARGUMENT, "+", "-", "*");

-- exponentiation and De Moivre's rule :

  function "**" ( c : double_complex; m : integer ) return double_complex;

  function de_Moivre (n,j : natural; c : double_complex) return double_complex;

   -- DESCRIPTION :
   --   returns the j-th solution of the equation x^n-c=0

private

  ZERO: constant double_float := double_float(0);
  ONE:  constant double_float := double_float(1);

  type double_complex is
    record
        RE : double_float := ZERO;
        IM : double_float := ZERO;
    end record;

  I: constant double_complex := ( ZERO, ONE );

end Complex_Numbers;
