with text_io;                          use text_io;
with Complex_Numbers;                  use Complex_Numbers;
with Float_Vectors;
with Float_Vectors_of_Vectors;
with Solutions;                        use Solutions;

package Drivers_for_Path_Directions is

-- DESCRIPTION :
--   This package provides driver routines for the computation of the
--   direction of the solution paths.

  procedure Init_Path_Directions
               ( n,nv : in natural;
                 v : in out Float_Vectors_of_Vectors.Link_to_Vector;
                 errv : in out Float_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Initializes the data for path directions.

  -- ON ENTRY :
  --   n         dimension of the solution vectors;
  --   nv        number of solution paths.
 
  -- ON RETURN :
  --   v         nv zero vectors of range 1..n;
  --   errv      nv entries equal to 1.0.

  procedure Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj,report : in boolean;
                 v : in out Float_Vectors_of_Vectors.Vector;
                 errv : in out Float_Vectors.Vector;
                 target : in double_complex );

  -- DESCRIPTION :
  --   Performs the continuation with computation of path directions.

  procedure Write_Directions 
               ( file : in file_type; v : in Float_Vectors_of_Vectors.Vector;
                 errv : in Float_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the directions of the paths with their errors to the file.

end Drivers_for_Path_Directions;
