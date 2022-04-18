with Lists;
with Integer_Vectors;                use Integer_Vectors;
with Integer_Vectors_of_Vectors;
with Simplices;                      use Simplices;

package Triangulations is

-- DESCRIPTION :
--   This package exports a data structure for dealing with regular
--   triangulations of Newton polytopes.

-- DATA STRUCTURE :

  package Lists_of_Simplices is new Lists(Simplex);
  type Triangulation is new Lists_of_Simplices.List;

  Null_Triangulation : constant Triangulation 
                           := Triangulation(Lists_of_Simplices.Null_List);

-- CREATORS :

  function Create ( s : Simplex ) return Triangulation;

  -- DESCRIPTION :
  --   Returns a trianguation with one simplex.

  procedure Update ( t : in out Triangulation; x : in Link_to_Vector );
  procedure Update ( t : in out Triangulation; x : in Link_to_Vector;
                     newt : out Triangulation );

  -- DESCRIPTION :
  --   Computes new simplices which contain x, stores them in newt,
  --   and adds them to the triangulation.  

  procedure Update_One ( t : in out Triangulation; x : in Link_to_Vector );

  -- DESCRIPTION :
  --   Computes new simplices that contain x and adds them
  --   to the triangulation.
  --   Hereby it will be assumed that adding x causes no points of t to be
  --   interior.

  procedure Connect ( t : in out Triangulation );
  
  -- DESCRIPTION ;
  --   Given a collection of simplices, the appropiate connections
  --   between the simplices will be computed and made.

  procedure Flatten ( t : in out Triangulation );
   
  -- DESCRIPTION :
  --   All simplices in t will be flattened.

  -- REQUIRED :
  --   Cells that are already flattened are grouped at the end of the list.

-- SELECTORS :

  function Is_Vertex ( t : Triangulation; x : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the given vector x is a vertex of one of
  --   the simplices in t.

  function Vertices ( t : Triangulation )
                    return Integer_Vectors_of_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector with all vertices of the simplices in t.

  function Vertices ( t : Triangulation; x : Vector )
                    return Integer_Vectors_of_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector with all vertices of the simplices in t,
  --   which contain the given vector x. 

  function Is_In ( t : Triangulation; x : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the point x belongs to one of the simplices in t.

  function Is_In ( t : Triangulation; x : Vector ) return Simplex;

  -- DESCRIPTION :
  --   If the point belongs to one of the simplices in t, then this
  --   simplex is returned, otherwise, the Null_Simplex is returned.

  function Is_In ( t : Triangulation; s : Simplex ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the simplex s is contained in t.

  function Volume ( t : Triangulation ) return natural;

  -- DESCRIPTION :
  --   Computes n! times the volume of the simplices in the triangulation.

-- DESTRUCTOR :

  procedure Clear ( t : in out Triangulation );

  -- DESCRIPTION :
  --   Frees the allocated memory space.

end Triangulations;
