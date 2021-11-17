with Natural_Vectors;
with Complex_Numbers;                   use Complex_Numbers;

package body Laurent_to_Polynomial_Converters is

  function Laurent_Polynomial_to_Polynomial
             ( p : Complex_Multivariate_Laurent_Polynomials.Poly )
             return Complex_Multivariate_Polynomials.Poly is

    res : Complex_Multivariate_Polynomials.Poly;
    tt : Complex_Multivariate_Laurent_Polynomials.Term;

  begin
    Laurent_Polynomial_to_Polynomial(p,tt,res);
    Complex_Multivariate_Laurent_Polynomials.Clear(tt);
    return res;
  end Laurent_Polynomial_to_Polynomial;

  procedure Laurent_Polynomial_to_Polynomial
             ( l : in Complex_Multivariate_Laurent_Polynomials.Poly;
               t : out Complex_Multivariate_Laurent_Polynomials.Term; 
               p : out Complex_Multivariate_Polynomials.Poly ) is

    min : Complex_Multivariate_Laurent_Polynomials.Degrees 
           := Complex_Multivariate_Laurent_Polynomials.Minimal_Degrees(l);
    tt : Complex_Multivariate_Laurent_Polynomials.Term;

  begin
    for i in min'range loop
      min(i) := -min(i);
    end loop;
    tt.cf := CMPLX(1.0);
    tt.dg := min;
    p := Laurent_Polynomial_to_Polynomial(l,tt); t := tt;
  end Laurent_Polynomial_to_Polynomial;

  function Laurent_Polynomial_to_Polynomial
         ( l : Complex_Multivariate_Laurent_Polynomials.Poly;
           t : Complex_Multivariate_Laurent_Polynomials.Term )
         return Complex_Multivariate_Polynomials.Poly is

    res : Complex_Multivariate_Polynomials.Poly;
    use Complex_Multivariate_Laurent_Polynomials;

    procedure Laurent_Term_to_Term ( tt : in Term; cont : out boolean ) is

      rt : Complex_Multivariate_Polynomials.Term;

    begin
      rt.cf := tt.cf;
      rt.dg := new Natural_Vectors.Vector(tt.dg'range);
      for i in tt.dg'range loop
        rt.dg(i) := tt.dg(i) + t.dg(i);
      end loop;
      Complex_Multivariate_Polynomials.Plus_Term(res,rt);
      Complex_Multivariate_Polynomials.Clear(rt);
      cont := true;
    end Laurent_Term_to_Term;
    procedure LP2P is new Visiting_Iterator(Laurent_Term_to_Term);

  begin
    LP2P(l);
    return res;
  end Laurent_Polynomial_to_Polynomial;

  function Laurent_to_Polynomial_System ( p : Laur_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Laurent_Polynomial_to_Polynomial(p(i));
    end loop;
    return res;
  end Laurent_to_Polynomial_System;

end Laurent_to_Polynomial_Converters;
