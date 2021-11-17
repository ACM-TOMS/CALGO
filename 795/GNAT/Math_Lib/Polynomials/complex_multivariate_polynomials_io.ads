with text_io;                           use text_io;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package Complex_Multivariate_Polynomials_io is

-- DESCRIPTION :
--   This package contains routines for the input and output
--   of complex multivariate polynomials in symbolic form.

-- EXAMPLES : i denotes sqrt(-1)
--   ex1 : x**2*y + 1/2*z*y**2 - 2*z + y**3 + x - 1E9/-8.E-6* y + 3;
--   ex2 : x^2*y + z*y^2 - 2*z + y^3 + x - y + 3;
--   ex3 : (1.01 + 2.8*i)*x1**2*x2 + x3**2*x1 - 3*x1 + 2*x2*x3 - 3;
--   ex4 : (x1^2*x2 + x3^2*x1 - 3*x1 + 2*x2*x3 - 3)*x2**2*(x2-1+i);

-- DATASTRUCTURE AND CONSTANT NEEDED :

  type Power is ('*','^');
  delimiter : constant character := ';';

-- EXCEPTIONS :

  ILLEGAL_CHARACTER : exception;
      -- occurs when a character is found unexpected to be

  ILLEGAL_SYMBOL : exception;
      -- occurs when an illegal symbol is read

  ILLEGAL_OPERATION : exception;
      -- occurs when one tries to perform illegal operations with
      -- polynomials

  INFINITE_NUMBER : exception;
      -- occurs when a rational coefficient has a denominator = 0

  OVERFLOW_OF_UNKNOWNS : exception;
      -- occurs when the number of unknowns turns to be bigger
      -- than first was mentioned

  BAD_BRACKET : exception;
      -- occurs when a bracket, like '(' or ')' is misplaced

-- THE INPUT OPERATIONS :

  procedure get ( n : in natural; p : out Poly );
  procedure get ( file : in file_type; n : in natural; p : out Poly );
  procedure get ( p : out Poly );
  procedure get ( file : in file_type; p : out Poly );

  -- DESCRIPTION :
  --   A polynomial will be read from file or standard input.
  --   No ordening on the monomials is assumed.

  -- REQUIRED : 
  --  * all unknows must begin with a letter and may have
  --    no symbols like '+', '-', '*', '^', '/', ';' or brackets in them,
  --    i = sqrt(-1) is reserved for complex numbers representation;
  --  * each symbol is limited to 3 characters;
  --  * the input is terminated by the delimiter;
  --  * no blanks may occur in the numbers;
  --  * if specified, then the file must already be opened for input.

  -- NOTE :
  --   The end_of_line symbol is not read.

  -- ON ENTRY :
  --   file       file where the input is,
  --              if not specified, then standard input is assumed;
  --   n          the number of unknowns,
  --              if not specified, then n will first be read.

  -- ON RETURN :
  --   p          a polynomial in n unknowns.

-- THE OUTPUT OPERATIONS :

  procedure put ( n : out natural; p : in Poly; pow : in Power := '*' );
  procedure put ( file : in file_type; n : out natural; p : in Poly; 
                  pow : in Power := '*' );
  procedure put ( p : in Poly; pow : in Power );
  procedure put ( file : in file_type; p : in Poly; pow : in Power );
  procedure put ( p : in Poly );
  procedure put ( file : in file_type; p : in Poly );

  -- DESCRIPTION :
  --   A polynomial in n unknowns is written on file or on standard output.

  -- REQUIRED :
  --   If specified, then the file is must already been opened for output.

  -- ON ENTRY :
  --   file       file where the output must come,
  --              if not specified, then standard output is assumed;
  --   p          a polynomial in n unknows;
  --   pow        kind of power symbol used, the default is '*'.

  -- ON RETURN :
  --   n          the number of unknows of p,
  --              if not specified, n will first be written.

  procedure put_line ( p : in Poly );
  procedure put_line ( file : in file_type; p : in Poly );
  procedure put_line ( p : in Poly; pow : in Power );
  procedure put_line ( file : in file_type; p : in Poly; pow : in Power );

  -- DESCRIPTION :
  --   Every term of the polynomial p will be written on a separate line.

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays on screen the formatting rules as on-line help facility.

end Complex_Multivariate_Polynomials_io;
