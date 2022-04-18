with Integer_Vectors; 
with Complex_Numbers;                   use Complex_Numbers;

package body Polynomial_to_Laurent_Converters is

  function Polynomial_to_Laurent_Polynomial
             ( p : Complex_Multivariate_Polynomials.Poly )
             return Complex_Multivariate_Laurent_Polynomials.Poly is

    res : Complex_Multivariate_Laurent_Polynomials.Poly
        := Complex_Multivariate_Laurent_Polynomials.Null_Poly;

    use Complex_Multivariate_Polynomials;

    procedure Term_to_Laurent_Term ( t : in Term; cont : out boolean ) is

      rt : Complex_Multivariate_Laurent_Polynomials.Term;

    begin
      rt.cf := t.cf;
      rt.dg := new Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        rt.dg(i) := t.dg(i);
      end loop;
      Complex_Multivariate_Laurent_Polynomials.Plus_Term(res,rt);
      Complex_Multivariate_Laurent_Polynomials.Clear(rt);
      cont := true;
    end Term_to_Laurent_Term;
    procedure P2LP is new Visiting_Iterator(Term_to_Laurent_Term);

  begin
    P2LP(p);
    return res;
  end Polynomial_to_Laurent_Polynomial;

  function Polynomial_to_Laurent_System ( p : Poly_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Polynomial_to_Laurent_Polynomial(p(i));
    end loop;
    return res;
  end Polynomial_to_Laurent_System;

end Polynomial_to_Laurent_Converters;
