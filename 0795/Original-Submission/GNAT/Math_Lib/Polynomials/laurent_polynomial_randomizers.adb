with Random_Number_Generators;
with Complex_Numbers,Integer_Vectors;          use Complex_Numbers;

package body Laurent_Polynomial_Randomizers is

  function Real_Randomize ( p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Randomize_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := CMPLX(Random_Number_Generators.Random);
      rt.dg := new Integer_Vectors.Vector'(t.dg.all);
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
      rt.dg := new Integer_Vectors.Vector'(t.dg.all);
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
      rt.dg := new Integer_Vectors.Vector'(t.dg.all);
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
      rt.dg := new Integer_Vectors.Vector'(t.dg.all);
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Perturb_Term;
    procedure Perturb_Terms is new Visiting_Iterator (Perturb_Term);
  begin
    Perturb_Terms(p);
    return res;
  end Complex_Perturb;

  function Real_Randomize ( p : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Real_Randomize(p(i));
    end loop;
    return res;
  end Real_Randomize;

  function Complex_Randomize ( p : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Complex_Randomize(p(i));
    end loop;
    return res;
  end Complex_Randomize;

  function Complex_Randomize1 ( p : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Complex_Randomize1(p(i));
    end loop;
    return res;
  end Complex_Randomize1;

  function Complex_Perturb ( p : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(p'range);
  begin
    for i in res'range loop
      res(i) := Complex_Perturb(p(i));
    end loop;
    return res;
  end Complex_Perturb;

end Laurent_Polynomial_Randomizers;
