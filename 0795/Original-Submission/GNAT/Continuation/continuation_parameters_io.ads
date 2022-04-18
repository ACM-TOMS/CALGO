with text_io;                   use text_io;

package Continuation_Parameters_io is

-- DESCRIPTION :
--   This package provides the primitives for a user interface to determine
--   the numerical parameters in the polynomial continuation.

  procedure put;
  procedure put ( file : in file_type );

  -- DESCRIPTION :
  --   Produces an overview of all values of the continuation parameters.

  procedure get ( k : out natural );

  -- DESCRIPTION :
  --   Prompts for a number k of a continuation parameter to change.
  --   If k < 34, then the corresponding parameter will be changed.

end Continuation_Parameters_io;
