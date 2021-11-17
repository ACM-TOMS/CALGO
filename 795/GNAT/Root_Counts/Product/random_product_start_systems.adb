with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Random_Number_Generators;           use Random_Number_Generators;
with Natural_Vectors;                    use Natural_Vectors;
with Complex_Vectors,Complex_Numbers;    use Complex_Numbers;
with Integer_Vectors,Complex_Matrices;   use Complex_Matrices;
with Complex_Linear_System_Solvers;      use Complex_Linear_System_Solvers;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;

with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;        use Degrees_in_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Degree_Structure;                   use Degree_Structure;
with Random_Product_System;

package body Random_Product_Start_Systems is

  procedure Add_Hyperplanes ( i : in natural; p : in Poly;
                              z : in partition; m : in natural;
                              dg : in Natural_Vectors.Vector) is
  -- DESCRIPTION : 
  --   This routine adds hyperplanes to the Random Product System 
  --   according to the partition and the degrees of the polynomial 

  -- ON ENTRY :
  --   i        the number of the polynomial in the system;
  --   p        the i-the polynomial of the system;
  --   z        a partition of the set of unknowns of p;
  --   m        the number of sets in z;
  --   dg       the degrees of the polynomial in the sets of z,
  --            for j in 1..m : dg(j) = Degree(p,z(j)).

    n : natural := Number_Of_Unknowns(p);
    h : Complex_Vectors.Vector(0..n);

  begin
    for k in 1..m loop
      for l in 1..dg(k) loop
        h(0) := random1;
        for j in 1..n loop
          if Is_In(z(k),j)
           then h(j) := Random1;
           else h(j) := CMPLX(0.0);
          end if;
        end loop;
        Random_Product_System.Add_Hyperplane(i,h);
      end loop;
    end loop;
  end Add_Hyperplanes;

  procedure Construct_Random_Product_System ( p : in Poly_Sys ) is

  -- DESCRIPTION : 
  --   This procedure constructs a random product system by
  --   finding a good partition for each equation of the system p.

    n : natural := p'length;
    m : natural;
  begin
    if Degree_Structure.Empty
     then Degree_Structure.Create(p);
    end if;
    for i in p'range loop
      m := Degree_Structure.Get(i);
      declare
        z : partition(1..m);
        dg : Natural_Vectors.Vector(1..m);
      begin
        Degree_Structure.Get(i,z,dg);
        Add_Hyperplanes(i,p(i),z,m,dg);
        Clear(z);
      end;
    end loop;
  end Construct_Random_Product_System;

  procedure RPQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List; nl : out natural) is

    n : natural := p'length;

  begin
    Random_Product_System.Init(n);
    Construct_Random_Product_System(p);
    q := Random_Product_System.Polynomial_System;
    Random_Product_System.Solve(sols,nl);
    Degree_Structure.Clear;
    Random_Product_System.Clear;
  end RPQ;

  procedure Solve ( i,n : in natural; sols : in out Solution_List;
                    a : in out Matrix; b : in out Complex_Vectors.Vector;
                    acc : in out partition ) is
  begin
    if i > n
     then declare
            aa : Matrix(a'range(1),a'range(2));
            bb : Complex_Vectors.Vector(b'range);
            rcond : double_float;
            ipvt : Integer_Vectors.Vector(b'range);
          begin
            for k in a'range(1) loop
              for l in a'range(2) loop
                aa(k,l) := a(k,l);
              end loop;
              bb(k) := b(k);
            end loop;
            lufco(aa,n,ipvt,rcond);
           -- put("rcond : "); put(rcond); new_line;
            if rcond + 1.0 /= 1.0
             then lusolve(aa,n,ipvt,bb);
                  declare
                    s : Solution(n);
                  begin
                    s.t := CMPLX(0.0);
                    s.m := 1;
                    s.v := bb;
                    Add(sols,s);
                  end;
            end if;
          exception
            when numeric_error => return;
          end;
     else declare
            h : Complex_Vectors.Vector(0..n);
            count : natural := 0;
            z : partition(1..n);
            m : natural;
            d : Natural_Vectors.Vector(1..n);
          begin
            Degree_Structure.Get(i,z,d);
            m := Degree_Structure.Get(i);
            for j1 in 1..m loop
              if Degree_Structure.Admissible(acc,i-1,z(j1))
               then acc(i) := Create(z(j1));
                    for j2 in 1..d(j1) loop
                      count := count + 1;
                      h := Random_Product_System.Get_Hyperplane(i,count);
                      b(i) := -h(0);
                      for k in 1..n loop
                        a(i,k) := h(k);
                      end loop;
                      Solve(i+1,n,sols,a,b,acc);
                    end loop;
                    Clear(acc(i));
               else count := count + d(j1);
              end if;
            end loop;
            Clear(z);
          end;
    end if;
  end Solve;

  procedure Solve_Start_System 
                ( n : in natural; sols : in out Solution_List ) is

    m : Matrix(1..n,1..n);
    v : Complex_Vectors.Vector(1..n);
    acc : Partition(1..n);

  begin
    for i in 1..n loop
      for j in 1..n loop
        m(i,j) := CMPLX(0.0);
      end loop;
      v(i) := CMPLX(0.0);
    end loop;
    Solve(1,n,sols,m,v,acc);
  end Solve_Start_System;

  procedure GBQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List ) is

    n : natural := p'length;

  begin
    Random_Product_System.Init(n);
    Construct_Random_Product_System(p);
    q := Random_Product_System.Polynomial_System;
    Solve_Start_System(n,sols);
    Degree_Structure.Clear;
    Random_Product_System.Clear;
  end GBQ;

end Random_Product_Start_Systems;
