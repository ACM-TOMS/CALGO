with Complex_Numbers,Complex_Norms;        use Complex_Numbers,Complex_Norms;

procedure Durand_Kerner ( p : in Vector; z,res : in out Vector;
			  maxsteps : in natural; eps : in double_float;
			  nb : out natural ) is

  pp : Vector(p'range);

  function Horner ( p : Vector; x : double_complex ) return double_complex is

  -- DESCRIPTION :
  --   Returns (..((a[n]*x + a[n-1])*x + a[n-2])*x + .. + a[1])*x + a[0].

  begin
    if p'first = p'last
     then return p(p'first);
     else declare
	    res : double_complex := p(p'last);
          begin
	    for i in reverse p'first..(p'last-1) loop
              res := res*x + p(i);
            end loop;
	    return res;
          end;
    end if;
  end Horner;

  procedure DK ( p : in Vector; z,res : in out Vector ) is

  -- DESCRIPTION :
  --   Computes one step in the Durand-Kerner iteration.

    function Compute_q ( i : integer; a : Vector ) return double_complex is

    -- DESCRIPTION :
    --   Computes the quotient needed in the Durand-Kerner step.

      res : double_complex;

    begin
      res := CMPLX(1.0);
      for j in a'range loop
        if j /= i
         then res := res*(a(i)-a(j));
        end if;
      end loop;
      return res;
    end Compute_q;

  begin
    for i in z'range loop
      z(i) := z(i) - Horner(p,z(i))/Compute_q(i,z);
      res(i) := Horner(p,z(i));
    end loop;
  end DK;

begin
  if p(p'last) /= CMPLX(1.0)
   then for i in p'range loop
	  pp(i) := p(i) / p(p'last);
        end loop;
   else for i in p'range loop
	  pp(i) := p(i);
        end loop;
  end if;
  for k in 1..maxsteps loop
    nb := k;
    DK(pp,z,res);  
    Write(k,z,res);
    exit when (Norm2(res) <= eps);
  end loop;
end Durand_Kerner;
