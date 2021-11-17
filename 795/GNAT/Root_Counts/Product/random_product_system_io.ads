with text_io;    use text_io;

package Random_Product_System_io is

-- DESCRIPTION :
--   This package contains some routines for input and output
--   of random product systems.

  procedure get ( n : out natural );
  procedure get ( file : in file_type; n : out natural );

  -- DESCRIPTION :
  --   This procedure reads from the standard input or from file the
  --   information necessary to create a random product structure.

  -- ON ENTRY (standard input or file) :
  --   n        dimension of the problem
  --   polynomial to initialize the Symbol_Table
  --   the random product system:
  --     each equation consist of linear polynomials,
  --     the null polynomial ends the list for each equation.

  -- ON RETURN :
  --   n        dimension of the problem
  --   the internal data of the package RPS are filled in.

  procedure put ( n,fore,after,exp : in natural );
  procedure put ( file : in file_type; n,fore,after,exp : in natural );

  -- DESCRIPTION :
  --   This procedure writes the hyperplanes on the standard input or on file.

  -- ON ENTRY :
  --   n        dimension of the problem;
  --   fore     number of digits for the `.' while printing floats;
  --   after    number of digits after the `.';
  --   exp      number of digits of the exponent.

end Random_Product_System_io;
