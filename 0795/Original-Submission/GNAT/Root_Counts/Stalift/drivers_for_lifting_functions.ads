with text_io;                           use text_io;
with Integer_Vectors_of_Vectors;
with Float_Vectors_of_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Float_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Float_Vector_Lists;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;

package Drivers_for_Lifting_Functions is

  function Read_Integer_Lifting
              ( l : Lists_of_Integer_Vectors.List )
              return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION : interactive point-wise lifting of a list.
  --   Returns the list of lifted points.

  function Read_Float_Lifting
              ( l : Lists_of_Float_Vectors.List )
              return Lists_of_Float_Vectors.List;

  -- DESCRIPTION : interactive point-wise lifting of a list.
  --   Returns the list of lifted points.

  procedure Driver_for_Lifting_Functions
              ( file : in file_type; p : in Poly_Sys;
                points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lilifu : in out Integer_Vectors_of_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Displays the menu with the available lifting functions and
  --   performs the selected integer lifting function.

  -- NOTE : it is assumed that different supports are submitted.

  -- ON ENTRY :
  --   file     file that must be opened for output;
  --   p        polynomial system;
  --   points   supports of the system p.

  -- ON RETURN :
  --   lifted   the lifted support sets;
  --   lilifu   vectors used for linear lifting, otherwise lilifu = null.

  procedure Driver_for_Lifting_Functions
              ( file : in file_type; p : in Poly_Sys;
                points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
                lilifu : in out Float_Vectors_of_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Displays the menu with the available lifting functions and
  --   performs the selected floating-point lifting function.

  -- NOTE : it is assumed that different supports are submitted.

  -- ON ENTRY :
  --   file     file that must be opened for output;
  --   p        polynomial system;
  --   points   supports of the system p.

  -- ON RETURN :
  --   lifted   the lifted support sets;
  --   lilifu   vectors used for linear lifting, otherwise lilifu = null.

  procedure Driver_for_Lifting_Functions
              ( file : in file_type; p : in Poly_Sys;
                ipoints : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                fltlif : out boolean;
                fpoints : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
                ilifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
                ililifu : in out Integer_Vectors_of_Vectors.Link_to_Vector;
                flilifu : in out Float_Vectors_of_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   The user has the choice for integer or floating-point lifting.
  --   On return, output parameter fltlif is true if the user wants
  --   floating-point lifting and false otherwise.
  --   Depending on fltlif, the appropriate parameters are determined.

end Drivers_for_Lifting_Functions;
