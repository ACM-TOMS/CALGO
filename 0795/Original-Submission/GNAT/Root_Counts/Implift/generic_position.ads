with Complex_Polynomial_Systems;   use Complex_Polynomial_Systems;
with Trees_of_Vectors;             use Trees_of_Vectors;

function Generic_Position ( p : Poly_Sys; tv : Tree_of_Vectors )
                          return boolean;

-- DESCRIPTION :
--   Returns true if the polytopes of p are in generic position w.r.t. tv.
--   If false, then it is undecided, because tv contains not all candidates.

-- REQUIRED :
--   tv has been computed by recursive procedure for the mixed volume.

-- ON ENTRY :
--   p         polynomial system;
--   tv        tree of vectors with useful directions.
