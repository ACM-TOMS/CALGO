with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Complex_Numbers,Natural_Vectors;  use Complex_Numbers,Natural_Vectors;
with Random_Number_Generators;         use Random_Number_Generators;

package body Homogenization is

  function Homogeneous_Part ( p : Poly ) return Poly is

    res : Poly := Null_Poly;
    d : integer := Degree(p);

    procedure Homogeneous_Term ( t : in Term; continue : out boolean ) is
    begin
      if Sum(t.dg) = d
       then Plus_Term(res,t);
            continue := true;
       else continue := false;
      end if;
    end Homogeneous_Term;
    procedure Homogeneous_Terms is new Visiting_Iterator(Homogeneous_Term);

  begin
    Homogeneous_Terms(p);
    return res;
  end Homogeneous_Part;

  function Homogeneous_Part ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Homogeneous_Part(p(i));
    end loop;
    return res;
  end Homogeneous_Part;

  function Real_Random_Hyperplane ( n : natural ) return Poly is

    res : Poly;
    t : Term;

  begin
    t.dg := new Natural_Vectors.Vector'(1..n => 0);
    t.cf := CMPLX(Random);
    res := Create(t);
    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..n loop
      t.dg := new Natural_Vectors.Vector'(1..n => 0);
      t.dg(i) := 1;
      t.cf := CMPLX(Random);
      Plus_Term(res,t);
      Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Complex_Random_Hyperplane ( n : natural ) return Poly is

    res : Poly;
    t : Term;

  begin
    t.dg := new Natural_Vectors.Vector'(1..n => 0);
    t.cf := Random;
    res := Create(t);
    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..n loop
      t.dg := new Natural_Vectors.Vector'(1..n => 0);
      t.dg(i) := 1;
      t.cf := Random;
      Plus_Term(res,t);
      Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Standard_Hyperplane ( n,i : natural ) return Poly is

  -- DESCRIPTION : Returns x_i - 1.

    res : Poly;
    t : Term;
  begin
    t.dg := new Natural_Vectors.Vector'(1..n => 0);
    t.cf := CMPLX(-1.0);
    res := Create(t);
    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Natural_Vectors.Vector'(1..n => 0);
    t.dg(i) := 1;
    t.cf := CMPLX(1.0);
    Plus_Term(res,t);
    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  procedure Construct_Real_Random_Hyperplanes
                ( s : in out Poly_Sys; m : natural ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
                ( s : in out Poly_Sys; m : natural ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Standard_Hyperplanes ( s : in out Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : natural := s'length;

  begin
    for i in s'range loop
      Clear(s(i));
      s(i) := Standard_Hyperplane(n,i);
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Enlarge_Before ( p : in out Poly; m : in natural ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

   procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

     d : Degrees := new Natural_Vectors.Vector(1..(t.dg'last+m));

   begin
     for i in 1..m loop
       d(i) := 0;
     end loop;
     for i in t.dg'range loop
       d(i+m) := t.dg(i);
     end loop;
     Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
     t.dg := d;
     continue := true;
   end Enlarge_Term;
   procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_After ( p : in out Poly; m : in natural ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

   procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

     d : Degrees := new Natural_Vectors.Vector(1..(t.dg'last+m));

   begin
     for i in t.dg'range loop
       d(i) := t.dg(i);
     end loop;
     for i in (t.dg'last+1)..m loop
       d(i) := 0;
     end loop;
     Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
     t.dg := d;
     continue := true;
   end Enlarge_Term;
   procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  function Add_Equations ( s1 : Poly_Sys; s2 : Poly_Sys ) return Poly_Sys is

    n1 : natural := s1'last - s1'first + 1;
    n2 : natural := s2'last - s2'first + 1;
    res : Poly_Sys(1..(n1+n2));
    m : natural;

  begin
    for i in 1..n1 loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m;
            Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..n2 loop
      Copy(s2(i),res(n1+i));
      m := Number_Of_Unknowns(res(n1+i));
      if m < (n1+n2)
       then m := n1 + n2 - m;
            Enlarge_Before(res(n1+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equation ( s : Poly_Sys; p : Poly ) return Poly_Sys is

    n : natural := s'last - s'first + 1;
    m : natural;
    res : Poly_Sys(1..(n+1));

  begin
    for i in 1..n loop
      Copy(s(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < n+1
       then m := n + 1 - m;
            Enlarge_After(res(i),m);
      end if;
    end loop;
    m := Number_Of_Unknowns(p);
    if m < (n+1)
     then m := n + 1 - m;
          Enlarge_Before(res(n+1),m);
    end if;
    return res;
  end Add_Equation;

  function Add_Random_Hyperplanes
               ( s : Poly_Sys; m : natural; re : boolean ) return Poly_Sys is

    s2 : Poly_Sys(1..m);
    n : natural := s'last - s'first + 1;
    res : Poly_Sys(1..(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Standard_Hyperplanes
               ( s : Poly_Sys; m : natural ) return Poly_Sys is

    n : natural := s'length;
    res : Poly_Sys(1..(n+m));
    s2 : Poly_Sys(1..m);

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

end Homogenization;
