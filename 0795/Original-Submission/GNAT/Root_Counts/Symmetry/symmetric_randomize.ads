with Symmetry_Group; use Symmetry_Group;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;

function Symmetric_Randomize ( p : Laur_Sys; v,w : List_of_Permutations )
                             return Laur_Sys;

-- DESCRIPTION :
--   This function returns a polynomial system with the same Newton
--   polytopes, but with randomly generated coefficients, w.r.t. to
--   the symmetry relations given by v and w.

-- ON ENTRY :
--   p         a polynomial system with n equations in n unknowns;
--   v,w       representations of the symmetry group.

-- ON RETURN :
--   a symmetric randomized polynomial system.
