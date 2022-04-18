with Floating_Point_Numbers;     use Floating_Point_Numbers;
with Complex_Numbers;            use Complex_Numbers;
with Natural_Vectors;

package body Projective_Transformations is

  function Projective_Transformation ( p : Poly ) return Poly is
  
    deg : integer := Degree(p);
    res : Poly := Null_Poly;

    procedure Homogeneous_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural := 0;

    begin
      ht.cf := t.cf;
      ht.dg := new Natural_Vectors.Vector(t.dg'first..t.dg'last+1);
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        ht.dg(i) := t.dg(i);
      end loop;
      ht.dg(ht.dg'last) := deg-sum;
      Plus_Term(res,ht);
      Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(ht.dg));
      continue := true;
    end Homogeneous_Term;
    procedure Homogeneous_Terms is new Visiting_Iterator(Homogeneous_Term);

  begin
    Homogeneous_Terms(p);
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( p : in out Poly ) is
  
    res : Poly := Projective_Transformation(p);

  begin
    Copy(res,p); Clear(res);
  end Projective_Transformation;

  function Projective_Transformation ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation ( p : in out Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  function Projective_Transformation ( s : Solution ) return Solution is

    n : natural := s.n;
    r : Solution(n+1);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := CMPLX(1.0);
    r.t := s.t;
    r.m := s.m;
    return r;
  end Projective_Transformation;

  procedure Projective_Transformation ( sols : in out Solution_List ) is
  begin
    if Is_Null(sols)
     then null;
     else declare
            temp : Solution_List := sols;
            n : natural := Head_Of(sols).n;
            l : Link_To_Solution;
            s : Solution(n);
            s2 : Solution(n+1);
          begin
            while not Is_Null(temp) loop
              l := Head_Of(temp);
              s := l.all;
              s2 := Projective_Transformation(s);
              Clear(l);
              l := new Solution'(s2);
              Set_Head(temp,l);
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
  end Projective_Transformation;

  function Affine_Transformation ( s : Solution ) return Solution is

    n : natural := s.n;
    r : Solution(n-1);

  begin
    for i in 1..(n-1) loop
      if modulus(s.v(n)) + CMPLX(1.0) = CMPLX(1.0)
       then r.v(i) := CMPLX(10.0**10);
       else r.v(i) := s.v(i) / s.v(n);
      end if;
     end loop;
     r.t := s.t;
     r.m := s.m;
     return r;
  exception
    when numeric_error => r.v(1..(n-1)) := (1..(n-1) => CMPLX(10.0**10));
                          return r;
  end Affine_Transformation;

  procedure Affine_Transformation ( sols : in out Solution_List ) is
  begin
    if Is_Null(sols)
     then null;
     else declare
            n : natural := Head_Of(sols).n;
            s1 : Solution(n);
            s2 : Solution(n-1);
            temp : Solution_List := sols;
            l : Link_To_Solution;
          begin
            while not Is_Null(temp) loop
              l := Head_Of(temp);
              s1 := l.all;
              s2 := Affine_Transformation(s1);
              Clear(l);
              l := new Solution'(s2);
              Set_Head(temp,l);
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
  end Affine_Transformation;

end Projective_Transformations;
