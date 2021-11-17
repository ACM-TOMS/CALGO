with Integer_Matrices;
with Complex_Numbers;                   use Complex_Numbers;
with Random_Number_Generators;          use Random_Number_Generators;
with Complex_Vectors,Complex_Matrices;  use Complex_Vectors,Complex_Matrices;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;       use Degrees_in_Sets_of_Unknowns;
with Random_Product_System;

package body m_Homogeneous_Start_Systems is

  procedure Create_Random_Hyperplanes ( index,n,d : in natural; s : in Set ) is
  begin
    for i in 1..d loop
      declare
        h : Complex_Vectors.Vector(0..n);
      begin
        h(0) := Random1;
        for j in 1..Dimension(s) loop
          if Is_In(s,j)
           then h(j) := random1;
           else h(j) := CMPLX(0.0);
          end if;
        end loop;
        Random_Product_System.Add_Hyperplane(index,h);
      end;
    end loop;
  end Create_Random_Hyperplanes;

  procedure Create_Random_System 
              ( n,m : natural; z : partition; d : Integer_Matrices.matrix ) is

  begin
    for j in 1..m loop
      for i in 1..n loop
        Create_Random_Hyperplanes(i,n,d(i,j),z(j));
      end loop;
    end loop;
  end Create_Random_System;

  procedure m_Homogeneous_Start_System
                 ( p : in Poly_Sys; z : in partition;
                   q : out Poly_Sys; qsols : in out Solution_List ) is

    n : constant natural := p'length;
    m : constant natural := z'last;
    d : constant Integer_Matrices.matrix := Degree_Table(p,z);
    nl : natural := 0;

  begin
    Random_Product_System.Init(n);
    Create_Random_System(n,m,z,d);
    Random_Product_System.Solve(qsols,nl);
    q := Random_Product_System.Polynomial_System;
  end m_Homogeneous_Start_System;

end m_Homogeneous_Start_Systems;
