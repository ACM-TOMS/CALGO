with Communications_with_User;           use Communications_with_User;
with Floating_Point_Numbers;             use Floating_Point_Numbers;
with integer_io,Numbers_io;              use integer_io,Numbers_io;
with Natural_Vectors;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers,Complex_Vectors;
with Complex_Numbers_io;                 use Complex_Numbers_io;
with Random_Number_Generators;           use Random_Number_Generators;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Interpolating_Homotopies;           use Interpolating_Homotopies;

procedure Driver_for_Interpolating_Homotopies
              ( file : in file_type; p : in Poly_Sys; z : in Partition;
                b : in out natural; q : out Poly_Sys; 
                qsols : in out Solution_List ) is

  n : constant natural := p'last;
  interpols : Solution_List;
  ib,scalind : natural;
  ans : character;

  function Random_Interpolating ( n,m : natural ) return Solution_List is

  -- DESCRIPTION :
  --   A list of m random n-dimensional vectors will be returned.
  --   The complex numbers will all have modulus one.

    res,res_last : Solution_List;
    s : Solution(n);

  begin
    s.t := CMPLX(0.0);
    s.m := 1;
    s.err := 0.0; s.rco := 1.0; s.res := 0.0;
    for i in 1..m loop
      for j in 1..n loop
        s.v(j) := random1;
      end loop;
      Append(res,res_last,s);
    end loop;
    return res;
  end Random_Interpolating;

  procedure Random_Linear_Scaler ( n : in natural; p : in Poly;
                                   v : out vector; l : out natural ) is

  -- DESCRIPTION :
  --   Returns a random vector of dimension n+1, with range 0..n.
  --   There will be a nonzero entry only for those unknowns that occur in p.

    res : Vector(0..n);
    last : natural := 0;

  begin
    for i in res'range loop
      if Degree(p,i) > 0
       then res(i) := random1;
            last := last + 1;
       else res(i) := CMPLX(0.0);
      end if;
    end loop;
    v := res; l := last;
  end Random_Linear_Scaler;

  function Scale ( sc,v : Vector; last : integer ) return double_complex is

  -- DESCRIPTION :
  --   Returns the last component of the vector v, that is v(last),
  --   such that sc(0) + sum sc(i)*v(i), i in v'range, holds.

    res : double_complex := sc(0);

  begin
    for i in v'first..last-1 loop
      res := res + sc(i)*v(i);
    end loop;
    res := -res/sc(last);
    return res;
  end Scale;

  function Random_Interpolating
                ( n,m : natural; scaler : Vector; scallast : natural ) 
                return Solution_List is

  -- DESCRIPTION :
  --   A list of m random n-dimensional vectors will be returned.
  --   The complex numbers will all have modulus one, except the last one,
  --   indicated by scallast, that has been chosen to satisfy the scaler 
  --   equation, defined by sum of scaler(i)*x_i = 0, with x_0 = 1.

    res,res_last : Solution_List;
    s : Solution(n);

  begin
    s.t := CMPLX(0.0);
    s.m := 1;
    for i in 1..m loop
      for j in 1..(n-1) loop
        s.v(j) := random1;
      end loop;
      s.v(n) := Scale(scaler,s.v,scallast);
      Append(res,res_last,s);
    end loop;
    return res;
  end Random_Interpolating;

  function Create ( v : Vector ) return Poly is

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Natural_Vectors.Vector'(v'first+1..v'last => 0);
    for i in v'range loop
      t.cf := v(i);
      if i > v'first
       then t.dg(i) := 1;
      end if;
      Plus_Term(res,t);
      if i > v'first
       then t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Create;

  function Interpolating_by_User ( n,m : natural ) return Solution_List is

  -- DESCRIPTION :
  --   A list of m n-dimensional vectors will be read from standard input.

    res,res_last : Solution_List;
    s : Solution(n);
   -- f1,f2 : double_float;

  begin
    put("Reading "); put(m,1); put(" "); put(n,1);
    put_line("-dimensional complex vectors.");
    for i in 1..m loop
      s.t := CMPLX(0.0);
      s.m := 1;
      s.err := 0.0; s.rco := 1.0; s.res := 0.0;
      put("Give the components of vector "); put(i,1); 
      put_line(" :");
      for j in 1..n loop
       -- Read_Double_Float(f1);
       -- Read_Double_Float(f2);
       -- s.v(j) := CMPLX(f1,f2);
        get(s.v(j));
      end loop;
      Append(res,res_last,s);
    end loop;
    return res;
  end Interpolating_by_User;

  procedure Driver_for_Interpolation is

  -- DESCRIPTION : interpolation without a scaling equation

    dp,ip : Poly_Sys(p'range);

  begin
    dp := Dense_Representation(p,z);
    ip := Independent_Representation(dp);
    ib := Independent_Roots(ip);
    if ib > b
     then ib := b;
    end if;
    put("The number of independent roots : "); put(ib,1); new_line;
    put("Do you want to give interpolation vectors by yourself ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then interpols := Interpolating_by_User(n,ib);
     else put_line("Random interpolating vectors will be generated.");
          interpols := Random_Interpolating(n,ib);
    end if;
    q := Interpolate(ip,ib,interpols);
    qsols := interpols; b := ib;
    Clear(dp); Clear(ip);
  end Driver_for_Interpolation;

  procedure Driver_for_Scaled_Interpolation is

  -- DESCRIPTION : interpolation with a scaling equation, p(scalind).

    dp,ip,pp,qq : Poly_Sys(p'range);
    scalvec : Vector(0..n);
    scalveclast : natural;

  begin
    Random_Linear_Scaler(n,p(scalind),scalvec,scalveclast);
    for i in p'range loop
      if i = scalind
       then pp(i) := Null_Poly;
       else pp(i) := p(i);
      end if;
    end loop;
    dp := Dense_Representation(pp,z);
    ip := Independent_Representation(dp);
    ib := Independent_Roots(ip,scalveclast);
    if ib > b
     then ib := b;
    end if;
    put("The number of independent roots : "); put(ib,1); new_line;
    put("Do you want to give interpolation vectors by yourself ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then interpols := Interpolating_by_User(n,ib);
     else put_line("Random interpolating vectors will be generated.");
          interpols := Random_Interpolating(n,ib,scalvec,scalveclast);
    end if;
    qq := Interpolate(ip,scalveclast,ib,interpols);
    qq(scalind) := Create(scalvec);
    qsols := interpols; b := ib; q := qq;
    Clear(dp); Clear(ip);
  end Driver_for_Scaled_Interpolation;

begin
  new_line;
  put("Give the number of the linear scaling equation (0 if none) : ");
  Read_Natural(scalind);
  if scalind = 0
   then Driver_for_Interpolation;
   else Driver_for_Scaled_Interpolation;
  end if;
end Driver_for_Interpolating_Homotopies;
