with Complex_Numbers;                    use Complex_Numbers;
with Random_Number_Generators;           use Random_Number_Generators;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Reduction_of_Polynomials;           use Reduction_of_Polynomials;

package body Reduction_of_Overconstrained_Systems is

  function Random_Square ( p : in Poly_Sys ) return Poly_Sys is

    m : constant natural := Number_of_Unknowns(p(p'first));
    res : Poly_Sys(1..m);

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in m+1..p'last loop
        declare
          a : double_complex := Random1;
          tmp : Poly := a*p(j);
        begin
          Plus_Poly(res(i),tmp);
          Clear(tmp);
        end;
      end loop;
    end loop;
    return res;
  end Random_Square;

  function Reduced_Square ( p : in Poly_Sys ) return Poly_Sys is

    m : constant natural := Number_of_Unknowns(p(p'first));
    res : Poly_Sys(1..m);

  begin
    for i in res'range loop
      Copy(p(i),res(i));
      for j in m+1..p'last loop
        declare
          tmp : Poly := Rpoly(res(i),p(j));
        begin
          Copy(tmp,res(i)); Clear(tmp);
        end;
      end loop;
    end loop;
    return res;
  end Reduced_Square;

end Reduction_of_Overconstrained_Systems;
