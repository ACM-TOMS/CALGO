with Random_Number_Generators;
with Complex_Numbers,Natural_Vectors;          use Complex_Numbers;
with Float_Matrices,Complex_Matrices;

package body Polynomial_Randomizers is

  function Real_Randomize ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := CMPLX(Random_Number_Generators.Random);
      rt.dg := new Natural_Vectors.Vector'(t.dg.all);
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator (Randomize_Term);
  begin
    Randomize_Terms(p);
    return res;
  end Real_Randomize;

  function Complex_Randomize ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := Random_Number_Generators.Random;
      rt.dg := new Natural_Vectors.Vector'(t.dg.all);
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator (Randomize_Term);
  begin
    Randomize_Terms(p);
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := Random_Number_Generators.Random1;
      rt.dg := new Natural_Vectors.Vector'(t.dg.all);
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Randomize_Term;
    procedure Randomize_Terms is new Visiting_Iterator (Randomize_Term);
  begin
    Randomize_Terms(p);
    return res;
  end Complex_Randomize1;

  function Complex_Perturb ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Perturb_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := t.cf + Random_Number_Generators.Random1;
      rt.dg := new Natural_Vectors.Vector'(t.dg.all);
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Perturb_Term;
    procedure Perturb_Terms is new Visiting_Iterator (Perturb_Term);
  begin
    Perturb_Terms(p);
    return res;
  end Complex_Perturb;

  function Real_Randomize ( p : Poly_Sys ) return Poly_Sys is
    res : Poly_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Real_Randomize(p(i));
    end loop;
    return res;
  end Real_Randomize;

  function Complex_Randomize ( p : Poly_Sys ) return Poly_Sys is
    res : Poly_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Complex_Randomize(p(i));
    end loop;
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Poly_Sys ) return Poly_Sys is
    res : Poly_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Complex_Randomize1(p(i));
    end loop;
    return res;
  end Complex_Randomize1;

  function Complex_Perturb ( p : Poly_Sys ) return Poly_Sys is
    res : Poly_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Complex_Perturb(p(i));
    end loop;
    return res;
  end Complex_Perturb;

  function Real_Random_Matrix ( n : natural )
                              return Float_Matrices.matrix is

  -- DESCRIPTION :
  --   Generates a dense n*n matrix with randomly generated floating
  --   point numbers.

    res : Float_Matrices.matrix(1..n,1..n);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Number_Generators.Random;
      end loop;
    end loop;
    return res;
  end Real_Random_Matrix;

  function Complex_Random_Matrix ( n : natural )
                                 return Complex_Matrices.matrix is

  -- DESCRIPTION :
  --   Generates a dense n*n matrix with randomly generated complex numbers.

    res : Complex_Matrices.matrix(1..n,1..n);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Number_Generators.Random;
      end loop;
    end loop;
    return res;
  end Complex_Random_Matrix;

  function Complex_Random_Matrix1 ( n : natural ) 
                                  return Complex_Matrices.matrix is

  -- DESCRIPTION :
  --   Generates a dense n*n matrix with randomly generated complex
  --   number with modulus 1.

    res : Complex_Matrices.matrix(1..n,1..n);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Number_Generators.Random1;
      end loop;
    end loop;
    return res;
  end Complex_Random_Matrix1;

  function Multiply ( m : Float_Matrices.matrix; 
                      p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Null_Poly;
      for j in p'range loop
        Plus_Poly(res(i),CMPLX(m(i,j))*p(j));
      end loop;
    end loop;
    return res;
  end Multiply;

  function Multiply ( m : Complex_Matrices.matrix;
                      p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in res'range loop
      res(i) := Null_Poly;
      for j in p'range loop
        Plus_Poly(res(i),m(i,j)*p(j));
      end loop;
    end loop;
    return res;
  end Multiply;

  function Randomize_with_Real_Matrix ( p : Poly_Sys ) return Poly_Sys is

    n : constant natural := p'length;
    mat : Float_Matrices.matrix(1..n,1..n) := Real_Random_Matrix(n);
    res : Poly_Sys(p'range) := Multiply(mat,p);

  begin
    return res;
  end Randomize_with_Real_Matrix;

  function Randomize_with_Complex_Matrix ( p : Poly_Sys ) return Poly_Sys is

    n : constant natural := p'length;
    mat : Complex_Matrices.matrix(1..n,1..n) := Complex_Random_Matrix(n);
    res : Poly_Sys(p'range) := Multiply(mat,p);

  begin
    return res;
  end Randomize_with_Complex_Matrix;

  function Randomize_with_Complex_Matrix1 ( p : Poly_Sys ) return Poly_Sys is

    n : constant natural := p'length;
    mat : Complex_Matrices.matrix(1..n,1..n) := Complex_Random_Matrix1(n);
    res : Poly_Sys(p'range) := Multiply(mat,p);

  begin
    return res;
  end Randomize_with_Complex_Matrix1;

end Polynomial_Randomizers;
