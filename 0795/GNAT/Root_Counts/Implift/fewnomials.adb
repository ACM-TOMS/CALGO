with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Integer_Vectors;
with Integer_Vectors_of_Vectors;

with Complex_Numbers;                   use Complex_Numbers;
with Complex_Vectors,Complex_Matrices;  use Complex_Vectors,Complex_Matrices;
with Complex_Linear_System_Solvers;     use Complex_Linear_System_Solvers;

with Complex_Multivariate_Laurent_Polynomials;
with Binomials;                         use Binomials;

package body Fewnomials is

  use Complex_Multivariate_Laurent_Polynomials;

  function Equal ( v : Integer_Vectors.Link_to_Vector; d : Degrees ) 
                 return boolean is
  begin
    return Integer_Vectors.Equal(v,Integer_Vectors.Link_to_Vector(d));
  end Equal;

  function Is_Fewnomial_System ( p : Laur_Sys ) return boolean is

    da : Integer_Vectors_of_Vectors.Vector(p'range);
    nb : natural;
    n : natural := p'last - p'first + 1;
    cnst : Integer_Vectors.Link_to_Vector;
    res : boolean;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
      found : boolean;
    begin
      found := false;
      for i in da'range loop
	exit when i > nb;
	if Equal(da(i),t.dg)
	 then found := true;
	      exit;
        end if;
      end loop;
      if not found
       then if nb < n
	     then nb := nb + 1;
		  found := true;
		  da(nb) := new Integer_Vectors.Vector(t.dg'range);
		  for j in t.dg'range loop
		    da(nb)(j) := t.dg(j);
                  end loop;
             elsif nb < n + 1
		 then nb := nb + 1;
		      cnst := new Integer_Vectors.Vector(t.dg'range);
		      found := true;
		      for j in t.dg'range loop
			cnst(j) := t.dg(j);
                      end loop;
		 elsif Equal(cnst,t.dg)
		     then found := true;
		     else res := false;
            end if;
      end if;
      cont := found;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator (Scan_Term);

  begin
    res := true;
    nb := 0;
    for i in p'range loop
      Scan_Terms(p(i));
      exit when not res;
    end loop;
    Integer_Vectors_of_Vectors.Clear(da);
    Integer_Vectors.Clear(cnst);
    return res;
  end Is_Fewnomial_System;

  procedure Scale ( cnst : in Integer_Vectors.Link_to_Vector;
                    da : in out Integer_Vectors_of_Vectors.Vector ) is

   -- DESCRIPTION :
   --   if cnst is not equal to the zero degrees,
   --   then all entries in da will be shifted along cnst.

    zero : boolean;

  begin
    zero := true;
    for i in cnst'range loop
      if cnst(i) /= 0
       then zero := false; exit;
      end if;
    end loop;
    if not zero
     then for i in da'range loop
	    Integer_Vectors.Min_Vector(da(i),cnst);
          end loop;
    end if;
  end Scale;

  procedure Initialize ( p : in Laur_Sys; n : in natural;
			 a : in out matrix; b : in out vector;
		         da : in out Integer_Vectors_of_Vectors.Vector;
                         fail : out boolean ) is

    -- DESCRIPTION :
    --   initializes the data needed for the solver;
    --   fail becomes true if the system p has too many different terms.

    cnst : Integer_Vectors.Link_to_Vector;
    da_numb,cnt : natural;
    fl : boolean;

    use Integer_Vectors;

    procedure Search_Term ( t : in Term; cont : out boolean ) is
      found : boolean;
    begin
      found := false;
      for i in da'range loop
	exit when i > da_numb;
	if Equal(da(i),t.dg)
	 then found := true;
	      a(cnt,i) := t.cf;
	      exit;
        end if;
      end loop;
      if not found
       then if da_numb < n
	     then da_numb := da_numb + 1;
		  da(da_numb) := new Integer_Vectors.Vector(t.dg'range);
		  for k in da(da_numb)'range loop
		    da(da_numb)(k) := t.dg(k);
                  end loop;
		  found := true;
		  a(cnt,da_numb) := t.cf;
	     elsif da_numb < n+1
		 then da_numb := da_numb + 1;
		      cnst := new Integer_Vectors.Vector(t.dg'range);
		      for k in cnst'range loop
			cnst(k) := t.dg(k);
                      end loop;
		      found := true;
		      b(cnt) := -t.cf;
                 elsif Equal(cnst,t.dg)
		     then found := true;
			  b(cnt) := -t.cf;
            end if;
      end if;
      if not found
       then cont := false;
	    fl := true;
       else cont := true;
	    fl := false;
      end if;
    end Search_Term;
    procedure Search is new Visiting_Iterator(Search_Term);

  begin
    da_numb := 0;
    for i in p'range loop
      for j in a'range(2) loop
	a(i,j) := CMPLX(0.0);
      end loop;
      b(i) := CMPLX(0.0);
      cnt := i;
      Search(p(i));
      exit when fl;
    end loop;
    if not fl and then (cnst /= null)
     then Scale(cnst,da);
	  Integer_Vectors.Clear(cnst);
	  fail := false;
     else fail := true;
    end if;
  end Initialize;

  procedure Solve ( p : in Laur_Sys; sols : in out Solution_List;
		    fail : out boolean ) is

    n : natural := p'last - p'first + 1;
    da : Integer_Vectors_of_Vectors.Vector(p'range);
    a : matrix(1..n,1..n);
    b : vector(1..n);
    fl : boolean;

  begin
    Initialize(p,n,a,b,da,fl);
    if fl
     then fail := true;
     else fail := false;
	  declare
	    piv : Integer_Vectors.Vector(1..n);
	    info : integer;
          begin
	    lufac(a,n,piv,info);
            if info = 0
	     then lusolve(a,n,piv,b);
		  fl := false;
		  for i in b'range loop
		    if modulus(b(i)) + 1.0 = 1.0
		     then fl := true;
			  exit;
                    end if;
                  end loop;
		  if not fl
                   then Binomials.Solve(da,b,n,sols);
                  end if;
            end if;
          end;
    end if;
    Integer_Vectors_of_Vectors.Clear(da);
  end Solve;

end Fewnomials;
