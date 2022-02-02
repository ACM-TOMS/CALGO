with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;
with Solutions;                       use Solutions;
with Partitions_of_Sets_of_Unknowns;  use Partitions_of_Sets_of_Unknowns;

package m_Homogeneous_Start_Systems is

-- DESCRIPTION :
--   the purpose of this package is to provide a routine
--   for the automatic construction of a m-homomogeneous
--   start system, given a partition of the set of unknowns.

  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in Partition;
                   q : out Poly_Sys; qsols : in out Solution_List );

  -- ON ENTRY :
  --   p           polynomial system;
  --   z           partition of the set of unknowns.

  -- ON RETURN :
  --   q           an m-homogeneous start system;
  --   qsols       the solutions of the start system q.

end m_Homogeneous_Start_Systems;
