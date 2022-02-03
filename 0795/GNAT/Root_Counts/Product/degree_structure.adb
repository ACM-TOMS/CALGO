with unchecked_deallocation,generate;
with text_io;                           use text_io;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers,Complex_Matrices;  use Complex_Numbers,Complex_Matrices;
with Complex_Numbers_io;                use Complex_Numbers_io;
with Complex_Linear_System_Solvers;     use Complex_Linear_System_Solvers;
with Degrees_in_Sets_of_Unknowns;       use Degrees_in_Sets_of_Unknowns;
with Random_Number_Generators;          use Random_Number_Generators;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Integer_Vectors;

package body Degree_Structure is

-- DECLARATIONS :

  type zd(m : natural) is record
    z : Partition(1..m);
    d : Vector(1..m);
  end record;
  type Link_To_zd is access zd;
  procedure free is new unchecked_deallocation(zd,Link_To_zd);

  type dgst is array(positive range <>) of Link_To_zd;
  type Link_To_dgst is access dgst;
  procedure free is new unchecked_deallocation(dgst,Link_To_dgst);

-- INTERNAL DATA :

  ds : Link_To_dgst;

-- CREATORS :

  procedure Find_Partition ( p : in Poly; 
                             z : in out Partition; m : in out natural; 
                             dg : in out Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   This routine finds an good partition for the polynomial p 

  -- ON ENTRY :
  --   p        a polynomial.

  -- ON RETURN :
  --   z        a partition of the set of unknowns of p;
  --   m        the number of sets in the partition z;
  --   dg       the degrees of the polynomial p in the sets of z,
  --            dg(i) = Degree(p,z(i)).

    n : natural := Number_Of_Unknowns(p);
    di : integer;
    added : boolean;

  begin
    for i in 1..n loop
      di := Degree(p,i);
      if di > 0
       then added := false;
            for j in 1..m loop
              if di = dg(j)
               then Add(z(j),i);
                    if Degree(p,z(j)) = dg(j)
                     then added := true;
                     else Remove(z(j),i);
                    end if;
              end if;
              exit when added;
            end loop;
            if not added
             then m := m + 1;
                  Add(z(m),i);
                  dg(m) := di;
            end if;
      end if;
    end loop;
  end Find_Partition;

  procedure Create ( p : in Poly_Sys ) is

    n : natural := p'length;
    z : Partition(1..n);
    d : Natural_Vectors.Vector(1..n) := (1..n => 0);
    m : natural;

  begin
    if ds /= null
     then Clear;
    end if;
    ds := new dgst(1..n);
    for i in p'range loop
      m := 0;
      Create(z,n);
      Find_Partition(p(i),z,m,d);
      ds(i) := new zd(m);
      for j in 1..m loop
        ds(i).z(j) := Create(z(j));
        ds(i).d(j) := d(j);
      end loop;
      Clear(z);
    end loop;
  end Create;

  procedure Put ( p : in Poly_Sys;
                  i,m : in natural; z : in Partition ) is

    n : natural := p'length;

  begin
    if ds = null
     then ds := new dgst(1..n);
    end if;
    ds(i) := new zd(m);
    for j in 1..m loop
      ds(i).z(j) := Create(z(j));
      ds(i).d(j) := Degree(p(i),z(j));
    end loop;
  end Put;

-- SELECTORS :

  function Empty return boolean is
  begin
    return (ds = null);
  end Empty;

  function Get ( i : natural ) return natural is
  begin
    return ds(i).m;
  end Get;

  procedure Get ( i : in natural; z : in out Partition;
                  d : out Natural_Vectors.Vector ) is
  begin
    for j in 1..ds(i).m loop
      z(j) := Create(ds(i).z(j));
      d(j) := ds(i).d(j);
    end loop;
  end Get;

-- COMPUTING THE GENERALIZED BEZOUT NUMBER :

  function Matrix_Criterion ( z : Partition ) return boolean is

  -- DESCRIPTION : 
  --   This is the matrix criterion for testing if
  --   a product is admissible or not.

    n : natural := z'last - z'first +1;
    mat : matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    eps : constant double_float := 10.0**(-10);
    r : double_complex;
    rcond : double_float;

  begin
    for i in 1..n loop
      r := CMPLX(double_float(i+1));
      for j in 1..n loop
        if Is_In(z(i),j)
         then mat(i,j) := r;
              r := r*r;
         else mat(i,j) := CMPLX(0.0);
        end if;
      end loop;
    end loop;
    lufco(mat,n,ipvt,rcond);
    return (abs(rcond) > eps);
  exception
    when others => return false;
  end Matrix_Criterion;
  
  function Admissible ( z : Partition; n : natural ) return boolean is

    temp : Partition(1..n);
    admis : boolean := true;

  begin
    temp(1) := Create(z(1));
    for i in 2..(n-1) loop
      temp(i) := Create(z(i));
      admis := Admissible(temp,i,z(i+1));
      exit when not admis;
    end loop;
    Clear(temp);
    return admis;
  end Admissible;

  function Admissible ( z : Partition; n : natural; s : Set )
                      return boolean is
  begin
    for k in 1..n loop
      if not Admissible(z,k,n,s)
       then return false;
      end if;
    end loop;
    return true;
  end Admissible;

  function Admissible ( z : Partition; k,n : natural; s : Set )
                      return boolean is

    type arr is array (integer range <>) of boolean;
    admis : boolean := true;

    procedure check (a : in arr; continue : out boolean) is

      u : Set := Create(s);

    begin
      for i in a'range loop
        if a(i)
         then Union(u,z(i));
        end if;
      end loop;
      admis := ( Extent_Of(u) >= k+1 );
      continue := admis;
      Clear(u);
    end check;

    procedure gen is new generate(arr,check);

  begin
    gen(k,1,n);
    return admis;
  end Admissible;

  procedure Compute ( i,n,sum : in natural; res : in out natural;
                      z : in out Partition ) is
  begin
    if i > n
     then res := res + sum;
     else -- Pick out a set and check if it is allowed :
          for j in 1..ds(i).m loop
            if ds(i).d(j) /= 0 and then Admissible(z,i-1,ds(i).z(j))
             then z(i) := Create(ds(i).z(j));
                  Compute(i+1,n,sum*ds(i).d(j),res,z);
                  Clear(z(i));
            end if;
          end loop;
    end if;
  end Compute;

  function Generalized_Bezout_Number return natural is

    res : natural := 0;
    n : natural := ds'length;
    z : Partition(1..n);

  begin
    Compute(1,n,1,res,z);
    return res;
  end Generalized_Bezout_Number;

  function Generalized_Bezout_Number ( p : in Poly_Sys ) return natural is
  begin
    Create(p);
    return Generalized_Bezout_Number;
  end Generalized_Bezout_Number;

-- DESTRUCTOR :

  procedure Clear is
  begin
    if ds /= null
     then for i in ds'range loop
            free(ds(i));
          end loop;
          free(ds);
    end if;
  end Clear;

end Degree_Structure;
