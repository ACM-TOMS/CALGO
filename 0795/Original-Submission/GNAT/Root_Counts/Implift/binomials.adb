with Complex_Numbers,Complex_Norms;  use Complex_Numbers,Complex_Norms;
with Integer_Vectors;
with Transforming_Solutions;         use Transforming_Solutions;

--with text_io,integer_io;             use text_io,integer_io;
--with Transformations_io;             use Transformations_io;
--with Integer_Vectors_of_Vectors_io;  use Integer_Vectors_of_Vectors_io;

package body Binomials is

-- AN INTERMEDIATE ROUTINE :  --

  procedure Factorize ( v : in out Vector; n : in natural;
                        t : in out Transfo ) is

  -- DESCRIPTION :
  --   This routines factorizes the binomial system defined by the
  --   degrees in the Vector v;

  -- ON ENTRY :
  --   v          defines the degrees of binomial system
  --   n          the number of unknowns to be eliminated

  -- ON RETURN :
  --   v          the factorized binomial system
  --   t          the transformations used to factorize

    tt : Transfo;
    pivot : integer;
    tmp : Integer_Vectors.Link_to_Vector;

  begin

    t := Id(v'last);

    for i in 1..n loop

     -- SEARCH FOR PIVOT :
      pivot := 0;
      for j in i..v'last loop
	if v(j)(i) /= 0
	 then pivot := j;
              --put("Pivot : "); put(pivot,1); new_line;
              exit;
        end if;
      end loop;

      if pivot /= 0
       then -- INTERCHANGE IF NECESSARY
	    if pivot /= i
	     then tmp := v(i);
		  v(i) := v(pivot);
		  v(pivot) := tmp;
                  --put_line("The pivoted matrix :"); put(v);
	    end if;

	    tt := Rotate(v(i),i);
            Apply(tt,v(i));
            for j in (i+1)..v'last loop
	      Apply(tt,v(j));
              v(j).all(i) := 0;
            end loop;

            --put("Step "); put(i,1); put_line(" :");
            --put_line("The transformation :"); put(tt);
            --put_line("The transformed matrix :"); put(v);

            Mult1(t,tt);
            Clear(tt);
      end if;
    end loop;
  end Factorize;

-- AUXILIARY ROUTINES :

  procedure Update ( sols1,sols2 : in out Solution_List ) is

    tmp : Solution_List := sols2;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      declare
        nls : Link_to_Solution := new Solution'(ls.all);
      begin
        Construct(nls,sols1);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(sols2);
  end Update;

  procedure Solve ( v : in Vector; cv : in Complex_Vectors.Vector;
                    i,n : in natural; sols : in out Solution_List;
                    acc : in out Complex_Vectors.Vector ) is
    t : Transfo;
    pivot,d : integer;
    a : double_complex;
    c : Complex_Vectors.Vector(cv'range) := cv;
    rv : Complex_Vectors.Vector(cv'range);
    tmp : Integer_Vectors.Link_to_Vector;
    nv : Vector(v'range);
    temp : array(v'range) of integer;
    newsols : Solution_List;
 
  begin
    for j in i..v'last loop
      nv(j) := new Integer_Vectors.Vector'(v(j).all);
    end loop;
   -- SEARCH FOR PIVOT :
    pivot := 0;
    for j in i..v'last loop
      if nv(j)(i) /= 0
       then pivot := j;
            exit;
      end if;
    end loop;
    if pivot /= 0
     then -- INTERCHANGE IF NECESSARY :
          if pivot /= i
           then tmp := nv(i);
         	nv(i) := nv(pivot);
		nv(pivot) := tmp;
		a := c(i);
		c(i) := c(pivot);
		c(pivot) := a;
	  end if;
	  t := Rotate(nv(i),i);
          Apply(t,nv(i));
          for j in (i+1)..nv'last loop
	    Apply(t,nv(j));
            temp(j) := nv(j)(i);
	    nv(j)(i) := 0;
          end loop;
         -- COMPUTE THE SOLUTIONS :
          d := nv(i)(i);
	  if d < 0
	   then d := -d;
		rv(i) := CMPLX(1.0)/c(i);
           else rv(i) := c(i);
          end if;
	  for j in 1..d loop
	    a := de_Moivre(d,j,rv(i));
	    acc(i) := a;
	    for k in (i+1)..nv'last loop
	      rv(k) := c(k)*a**(-temp(k));
            end loop;
            if i < n
	     then Solve(nv,rv,i+1,n,newsols,acc);
	     else -- i = n
		  declare
		    ls : Link_to_Solution;
		    s : Solution(n);
                  begin
		    s.t := CMPLX(0.0);
		    s.m := 1;
                    s.v := acc;
		    ls := new Solution'(s);
		    Construct(ls,newsols);
                  end;
            end if;
            Transform(t,newsols);
            Update(sols,newsols);
          end loop;
          Clear(t);
    end if;
    for j in i..nv'last loop
      Integer_Vectors.Clear(nv(j));
    end loop;
  end Solve;

-- THE SOLVER :

  procedure Solve ( v : in Vector; cv : in Complex_Vectors.Vector;
		    n : in natural; sols : in out Solution_List ) is

    wv : Vector(v'range);   -- workspace
    acc : Complex_Vectors.Vector(cv'range); -- accumulator

  begin
    for i in v'range loop
      wv(i) := new Integer_Vectors.Vector'(v(i).all);
    end loop;
    Solve(wv,cv,v'first,v'last,sols,acc);
    Clear(wv);
  end Solve;

-- COMPUTING THE RESIDUALS :

  procedure Residuals ( v : in Vector; cv : in Complex_Vectors.Vector;
			n : in natural; sols : in Solution_List;
			res : out Complex_Vectors.Vector ) is
    nb : natural;
    x : Complex_Vectors.Vector(cv'range);
    tmp : Solution_List;

    function Eval ( v : Vector; cv,x : Complex_Vectors.Vector )
                  return Complex_Vectors.Vector is

    -- DESCRIPTION :
    --   Computes the value of the vector x in the system defined by
    --   v and cv.

      y : Complex_Vectors.Vector(x'range);

    begin
      for i in x'range loop
        y(i) := CMPLX(1.0);
	for j in v(i)'range loop
	  y(i) := y(i)*x(j)**v(i)(j);
        end loop;
	y(i) := y(i) - cv(i);
      end loop;
      return y;
    end Eval;

  begin
    tmp := sols;
    nb := 1;
    while not Is_Null(tmp) loop
      x := Head_Of(tmp).v;
      res(nb) := CMPLX(Norm2(Eval(v,cv,x)));
      tmp := Tail_Of(tmp);
      nb := nb + 1;
    end loop;
  end Residuals;

end Binomials;
