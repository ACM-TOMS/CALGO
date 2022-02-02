with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers;                   use Complex_Numbers;
with Natural_Vectors,Integer_Vectors;
with Complex_Vectors,Complex_Matrices;  use Complex_Vectors,Complex_Matrices;
with Complex_Linear_System_Solvers;     use Complex_Linear_System_Solvers;
with Polynomial_Randomizers;            use Polynomial_Randomizers;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;       use Degrees_in_Sets_of_Unknowns;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package body Interpolating_Homotopies is

-- AUXILIARY OPERATIONS  :

  function Initial_Term ( p : Poly ) return Term is

  -- DESCRIPTION :
  --   Returns the initial term of p.

    res : Term;

    procedure Init_Term ( t : Term; continue : out boolean ) is
    begin
      res := t;
      continue := false;
    end Init_Term;
    procedure Init_Terms is new Visiting_Iterator(Init_Term);

  begin
    Init_Terms(p);
    return res;
  end Initial_Term;

  procedure Eliminate_Term ( p : in out Poly; dg : in Degrees ) is

  -- DESCRIPTION :
  --   Eliminates the monomial in p that has the given exponent vector.

    c : double_complex := Coeff(p,dg);

  begin
    if c /= CMPLX(0.0)
     then declare
            t : Term;
          begin
            t.cf := c; t.dg := dg;
            Min_Term(p,t);
          end;
    end if;
  end Eliminate_Term;

  function Admitted ( t : Term; i : natural; z : partition;
                      d : Integer_Matrices.matrix ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the given term does not violate the m-homogeneous
  --   structure, i.e., if Degree(t,z(j)) <= d(i,j), for j in z'range.

  begin
    for j in z'range loop
      if Degree(t,z(j)) > d(i,j)
       then return false; 
      end if;
    end loop;
    return true; 
  end Admitted;

  function Dense_Representation ( p : Poly; i : natural; z : Partition;
                                  d : Integer_Matrices.matrix ) return Poly is

  -- DESCRIPTION :
  --   Returns the dense representation of the given polynomial p, known as
  --   p(i) of some polynomial system, of the given m-homogeneous structure.
  --   The returned polynomial has all its coefficients equal to one.

    res : Poly := Null_Poly;
    maxdegs,accu : Natural_Vectors.Vector(d'range(1));

    procedure Generate_Monomials
                  ( k : in natural; max : in Natural_Vectors.Vector;
                    acc : in out Natural_Vectors.Vector) is

    -- DESCRIPTION :
    --   Generates all monomials with exponents in a box, with upper bounds
    --   given by the vector max.  Every monomial that respects the 
    --   m-homogeneous structure will be added to the result res.
    --   The current unknown is indicated by the parameter k.

      t : Term;

    begin
      for j in 0..max(k) loop
        acc(k) := j;
        t.cf := CMPLX(1.0);
        t.dg := new Natural_Vectors.Vector'(acc);
        if Admitted(t,i,z,d)
         then Plus_Term(res,t);
        end if;
        Clear(t);
        if k < max'last
         then Generate_Monomials(k+1,max,acc);
        end if;
      end loop;
    end Generate_Monomials;

  begin
    for j in maxdegs'range loop
      maxdegs(j) := Degree(p,j); accu(j) := 0;
    end loop;
    Generate_Monomials(1,maxdegs,accu);
    return res;
  end Dense_Representation;

  function Evaluate ( t : Term; x : Complex_Vectors.Vector )
                    return double_complex is

    res : double_complex := CMPLX(1.0);

  begin
    for i in x'range loop
      for k in 1..t.dg(i) loop
        res := res*x(i);
      end loop;
    end loop;
    return res;
  end Evaluate;

  procedure Interpolation_System
                ( p : in Poly; b : in natural; sols : in Solution_List;
                  mat : out Complex_Matrices.matrix;
                  rhs : out Complex_Vectors.vector ) is

  -- DESCRIPTION :
  --   Returns the matrix and right hand side vector of the linear system
  --   that expresses the interpolationg conditions.

    m : Complex_Matrices.matrix(1..b,1..b);
    v : Complex_Vectors.vector(1..b);
    cnt : natural := 0;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      tmp : Solution_List := sols;
      s : Link_to_Solution;

    begin
      for row in m'range(1) loop        -- fill in column indicated by cnt
        s := Head_Of(tmp);
        if cnt = 0
         then v(row) := -t.cf*Evaluate(t,s.v);
         elsif cnt <= b
             then m(row,cnt) := Evaluate(t,s.v);
             else v(row) := v(row) - t.cf*Evaluate(t,s.v);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      cnt := cnt + 1;
      continue := true;
    end Scan_Term;
    procedure Scan_Poly is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Poly(p);
    mat := m; rhs := v;
  end Interpolation_System;

  procedure Interpolation_System
                ( p : in Poly; i,b : in natural; sols : in Solution_List;
                  mat : out Complex_Matrices.matrix;
                  rhs : out Complex_Vectors.vector ) is

  -- DESCRIPTION :
  --   Returns the matrix and right hand side vector of the linear system
  --   that expresses the interpolationg conditions.  Monomials with degree
  --   one in x_i will not appear in the interpolation matrix.

    m : Complex_Matrices.matrix(1..b,1..b);
    v : Complex_Vectors.vector(1..b);
    cnt : natural := 0;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is

      tmp : Solution_List := sols;
      s : Link_to_Solution;

    begin
      for row in m'range(1) loop        -- fill in column indicated by cnt
        s := Head_Of(tmp);
        if cnt = 0
         then v(row) := -t.cf*Evaluate(t,s.v);
         else if cnt <= b and t.dg(i) /= 1
               then m(row,cnt) := Evaluate(t,s.v);
               else v(row) := v(row) - t.cf*Evaluate(t,s.v);
              end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if cnt = 0 or t.dg(i) /= 1
       then cnt := cnt + 1;
      end if;
      continue := true;
    end Scan_Term;
    procedure Scan_Poly is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Poly(p);
    mat := m; rhs := v;
  end Interpolation_System;

  procedure Construct_Interpolant ( p : in out Poly;
                                    cv : in Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   With the coefficients of its monomials, the interpolant will be
  --   constructed.

    cnt : natural := 0;

    procedure Scan_Term ( t : in out Term; continue : out boolean ) is
    begin
      if cnt > cv'last
       then continue := false;
       else if cnt >= cv'first
             then t.cf := cv(cnt);
            end if;
            cnt := cnt + 1;
            continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Poly is new Changing_Iterator(Scan_Term);

  begin
    Scan_Poly(p);
  end Construct_Interpolant;

  procedure Construct_Interpolant ( p : in out Poly; i : in natural;
                                    cv : in Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   With the coefficients of its monomials, the interpolant will be
  --   constructed.  Monomials that have degree one in x_i will be ignored.

    cnt : natural := 0;

    procedure Scan_Term ( t : in out Term; continue : out boolean ) is
    begin
      if cnt > cv'last
       then continue := false;
       else if cnt = 0 or t.dg(i) /= 1
             then if cnt >= cv'first
                   then t.cf := cv(cnt);
                  end if;
                  cnt := cnt + 1;
            end if;
            continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Poly is new Changing_Iterator(Scan_Term);

  begin
    Scan_Poly(p);
  end Construct_Interpolant;

  function Interpolate ( p : Poly; b : natural; sols : Solution_List )
                       return Poly is

  -- DESCRIPTION :
  --   Constructs the interpolating polynomial with the same monomial
  --   structure as the given polynomial.

    res : Poly := Complex_Randomize(p);
    mat : Complex_Matrices.Matrix(1..b,1..b);
    rhs : Complex_Vectors.Vector(1..b);
    ipvt : Integer_Vectors.Vector(1..b);
    info : integer;

  begin
    Interpolation_System(res,b,sols,mat,rhs);
    for i in ipvt'range loop
      ipvt(i) := i;
    end loop;
    lufac(mat,b,ipvt,info);
    if info = 0
     then lusolve(mat,b,ipvt,rhs);
    end if;
    Construct_Interpolant(res,rhs);
    return res;
  end Interpolate;

  function Interpolate ( p : Poly; i,b : natural; sols : Solution_List )
                       return Poly is

  -- DESCRIPTION :
  --   Constructs the interpolating polynomial with the same monomial
  --   structure as the given polynomial.  The monomials that have degree 
  --   one in x_i will not appear in the interpolation matrix.

    res : Poly := Complex_Randomize(p);
    mat : Complex_Matrices.Matrix(1..b,1..b);
    rhs : Complex_Vectors.Vector(1..b);
    ipvt : Integer_Vectors.Vector(1..b);
    info : integer;

  begin
    Interpolation_System(res,i,b,sols,mat,rhs);
    for j in ipvt'range loop
      ipvt(j) := j;
    end loop;
    lufac(mat,b,ipvt,info);
    if info = 0
     then lusolve(mat,b,ipvt,rhs);
    end if;
    Construct_Interpolant(res,i,rhs);
    return res;
  end Interpolate;

-- TARGET ROUTINES :

  function Dense_Representation
              ( p : Poly_Sys; z : partition ) return Poly_Sys is

    d : constant Integer_Matrices.matrix := Degree_Table(p,z);

  begin
    return Dense_Representation(p,z,d);
  end Dense_Representation;

  function Dense_Representation
              ( p : Poly_Sys; z : partition; d : Integer_Matrices.matrix )
              return Poly_Sys is

    res : Poly_Sys(d'range(1));

  begin
    for i in res'range loop
      if p(i) /= Null_Poly
       then res(i) := Dense_Representation(p(i),i,z,d);
      end if;
    end loop;
    return res;
  end Dense_Representation;

  function Independent_Representation ( p : Poly_Sys ) return Poly_Sys is
 
    res : Poly_Sys(p'range);
    it : Term;

  begin
    Copy(p,res);
    for i in res'range loop
      if p(i) /= Null_Poly
       then it := Initial_Term(res(i));
            for j in res'range loop
              if j /= i and then (p(j) /= Null_Poly)
               then Eliminate_Term(res(j),it.dg);
              end if;
            end loop;
      end if;
    end loop;
    return res;
  end Independent_Representation;

  function Independent_Roots ( p : Poly_Sys ) return natural is

    ntp : natural := 0;
    min : natural := ntp;

  begin
    for i in p'first..p'last loop
      if p(i) /= Null_Poly
       then ntp := Number_of_Terms(p(i));
            if min = 0
             then min := ntp;
             elsif ntp < min
                 then min := ntp;
            end if;
      end if;
    end loop;
    if min = 0
     then return 0;
     else return (min-1);
    end if;
  end Independent_Roots;

  function Number_of_Terms ( p : Poly; i : natural ) return natural is

  -- DESCRIPTION :
  --   Returns the number of monomials of p, without those monomials that
  --   have degree one in x_i.

    cnt : natural := 0;

    procedure Count_Term ( t : in Term; continue : out boolean ) is
    begin
      if t.dg(i) /= 1
       then cnt := cnt+1;
      end if;
      continue := true;
    end Count_Term;
    procedure Count_Terms is new Visiting_Iterator(Count_Term);

  begin
    Count_Terms(p);
    return cnt;
  end Number_of_Terms;

  function Independent_Roots ( p : Poly_Sys; i : natural ) return natural is

    ntp : natural := 0;
    min : natural := ntp;

  begin
    for j in p'first..p'last loop
      if p(j) /= Null_Poly
       then ntp := Number_of_Terms(p(j),i);
            if min = 0
             then min := ntp;
             elsif ntp < min
                 then min := ntp;
            end if;
      end if;
    end loop;
    if min = 0
     then return 0;
     else return (min-1);
    end if;
  end Independent_Roots;

  function Interpolate ( p : Poly_Sys; b : natural; sols : Solution_List )
                       return Poly_Sys is

    res : Poly_Sys(p'range);

  begin 
    for i in p'range loop
      if p(i) /= Null_Poly
       then res(i) := Interpolate(p(i),b,sols);
      end if;
    end loop;
    return res;
  end Interpolate;

  function Interpolate ( p : Poly_Sys; i,b : natural; sols : Solution_List )
                       return Poly_Sys is

    res : Poly_Sys(p'range);

  begin 
    for j in p'range loop
      if p(j) /= Null_Poly
       then res(j) := Interpolate(p(j),i,b,sols);
      end if;
    end loop;
    return res;
  end Interpolate;

end Interpolating_Homotopies;
