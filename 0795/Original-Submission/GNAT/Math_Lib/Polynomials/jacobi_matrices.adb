with Complex_Numbers;          use Complex_Numbers;

package body Jacobi_Matrices is

-- CREATORS :

  function Create ( p : Poly_Sys ) return Jacobi is

    res : Jacobi(p'range,1..Number_of_Unknowns(p(p'first)));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Diff(p(i),j);
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( j : Jacobi ) return Eval_Jacobi is

    res : Eval_Jacobi(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        res(k,l) := Create(j(k,l));
      end loop;
    end loop;
    return res;
  end Create;

  procedure Create ( p : Poly_Sys;
                     j : out Eval_Coeff_Jacobi; m : out Mult_Factors ) is
  
    nb : constant natural := Number_of_Unknowns(p(p'first));
    nbk : natural;

  begin
    for k in p'range loop
      nbk := Number_of_Terms(p(k));
      for l in 1..nb loop
        declare
          mkl : Vector(1..nbk);
        begin
          Diff(p(k),l,j(k,l),mkl);
          m(k,l) := new Complex_Vectors.Vector'(mkl);
        end;
      end loop;
    end loop;
  end Create;

-- OPERATIONS :

  function Equal ( j1,j2 : Jacobi ) return boolean is
  begin
    for k in j1'range loop
      for l in j2'range loop
        if not Equal(j1(k,l),j2(k,l))
         then return false;
        end if;
      end loop;
    end loop;
    return true;
  exception
    when CONSTRAINT_ERROR => return false;
  end Equal;

  procedure Copy ( j1 : in Jacobi; j2 : in out Jacobi ) is
  begin
    Clear(j2);
    for k in j1'range(1) loop
      for l in j2'range(2) loop
        Copy(j1(k,l),j2(k,l));
      end loop;
    end loop;
  end Copy;

  function "+" ( j1,j2 : Jacobi ) return Jacobi is 

    res : Jacobi(j1'range(1),j1'range(2));

  begin
    Copy(j1,res);
    Plus_Jacobi(res,j2);
    return res;
  end "+";

  function "-" ( j : Jacobi ) return Jacobi is

    res : Jacobi(j'range(1),j'range(2));

  begin
    Copy(j,res);
    Min_Jacobi(res);
    return res;
  end "-";

  function "-" ( j1,j2 : Jacobi ) return Jacobi is

    res : Jacobi(j1'range(1),j1'range(2));

  begin
    Copy(j1,res);
    Min_Jacobi(res,j2);
    return res;
  end "-";

  procedure Plus_Jacobi ( j1 : in out Jacobi; j2 : in Jacobi ) is
  begin
    for k in j1'range(1) loop
      for l in j1'range(1) loop
        Plus_Poly(j1(k,l),j2(k,l));
      end loop;
    end loop;
  end Plus_Jacobi;

  procedure Min_Jacobi ( j  : in out Jacobi ) is
  begin
    for k in j'range(1) loop
      for l in j'range(1) loop
        Min_Poly(j(k,l));
      end loop;
    end loop;
  end Min_Jacobi;

  procedure Min_Jacobi ( j1 : in out Jacobi; j2 : in Jacobi ) is 
  begin
    for k in j1'range(1) loop
      for l in j1'range(1) loop
        Min_Poly(j1(k,l),j2(k,l));
      end loop;
    end loop;
  end Min_Jacobi;

  function Eval ( j : Jacobi; x : Vector ) return matrix is

    m : matrix(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        m(k,l) := Eval(Poly(j(k,l)),x);
      end loop;
    end loop;
    return m;
  end Eval;

  function Eval ( j : Eval_Jacobi; x : Vector ) return matrix is

    m : matrix(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        m(k,l) := Eval(Eval_Poly(j(k,l)),x);
      end loop;
    end loop;
    return m;
  end Eval;

  function Eval ( j : Eval_Coeff_Jacobi; m : Mult_Factors;
                  c : Complex_Vectors_of_Vectors.Vector; x : Vector )
                return Matrix is
 
    res : matrix(j'range(1),j'range(2));

  begin
    for k in j'range(1) loop
      declare
        cm : Vector(c(k)'range);
      begin
        for l in j'range(2) loop
          for i in cm'range loop
            cm(i) := m(k,l)(i)*c(k)(i);
          end loop;
          res(k,l) := Eval(Eval_Coeff_Poly(j(k,l)),cm,x);
        end loop;
      end;
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( j : in out Jacobi ) is
  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        Clear(j(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( j : in out Eval_Jacobi ) is
  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        Clear(j(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( j : in out Eval_Coeff_Jacobi ) is
  begin
    for k in j'range(1) loop
      for l in j'range(2) loop
        Clear(j(k,l));
      end loop;
    end loop;
  end Clear;

  procedure Clear ( m : in out Mult_Factors ) is
  begin
    for k in m'range(1) loop
      for l in m'range(2) loop
        Clear(m(k,l));
      end loop;
    end loop;
  end Clear;
  
end Jacobi_Matrices;
