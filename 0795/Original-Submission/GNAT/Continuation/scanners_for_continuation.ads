with text_io;           use text_io;
with Complex_Numbers;   use Complex_Numbers;

package Scanners_for_Continuation is

-- DESCRIPTION :
--   This package provides scanning procedures for reading the
--   input data for polynomial continuation from file.

  procedure Scan_Homotopy_Parameters
              ( file : in file_type; k : out natural; a : out double_complex );

  -- DESCRIPTION :
  --   Scans the file for the homotopy parameters k and a.
  --   If they are not found, then both output parameters are zero.

  procedure Scan_Continuation_Parameters ( file : in file_type );

  -- DESCRIPTION :
  --   Scans the file for the continuation parameters.
  --   If no continuation parameters are given, then nothing happens.

  procedure Scan_Output_Parameter ( file : in file_type; op : out natural );

  -- DESCRIPTION :
  --   Scans the file for the output parameter (op) of the continuation.
  --   If no output code is given, then nothing happens.

end Scanners_for_Continuation;
