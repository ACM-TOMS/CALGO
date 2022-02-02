with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;
with Solutions;                          use Solutions;

package Fewnomials is

-- DESCRIPTION :
--   This package contains a solver for Laurent polynomial systems,
--   where there are at most n+1 different terms.

  function Is_Fewnomial_System ( p : Laur_Sys ) return boolean;

  -- DESCRIPTION :
  --   returns true if p has at most n+1 different terms.

  procedure Solve ( p : in Laur_Sys; sols : in out Solution_List;
		    fail : out boolean );

  -- DESCRIPTION :
  --   tries to solve the fewnomial system p.

  -- ON ENTRY :
  --   p          a Laurent polynomial system.

  -- ON RETURN :
  --   sols       solutions of p, all solutions have only nonzero components;
  --   fail       true if p has more than n+1 different terms
  --              or if p has an infinite number of solutions.

end Fewnomials;
