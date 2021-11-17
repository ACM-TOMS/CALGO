with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers;                   use Complex_Numbers;
with Natural_Vectors,Complex_Matrices;  use Complex_Matrices;
with Complex_Linear_System_Solvers;     use Complex_Linear_System_Solvers;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Reduction_of_Polynomials;          use Reduction_of_Polynomials;

package body Reduction_of_Polynomial_Systems is

-- AUXALIARY DATA FOR LINEAR REDUCTION :

  mach_eps : constant double_float := 10.0**(-13);

  type Degrees_Array is array(positive range <>) of Degrees;
  type Terms_Array   is array(positive range <>) of Term;
  type Boolean_Array is array(positive range <>) of boolean;

-- AUXILIARY PROCEDURES FOR LINEAR REDUCTION :

  procedure Pop_First_Term ( p : in out Poly; t : in out Term ) is

  -- DESCRIPTION :
  --   The term on return is the leading term of p.
  --   This term is removed from p.

    procedure First_Term ( tt : in Term; continue : out boolean ) is
    begin
      Copy(tt,t);
      continue := false;
    end First_Term;
    procedure Get_First_Term is new Visiting_Iterator(First_Term);

  begin
    t.cf := CMPLX(0.0);
    Get_First_Term(p);
    if t.cf /= CMPLX(0.0)
     then Min_Term(p,t);
    end if;
  end Pop_First_Term;

  procedure Leading_Terms ( p : in out Poly_Sys; ta : in out Terms_Array ) is

  -- DESCRIPTION :
  --   Puts the leading terms of the polynomials in p in the array ta.
  --   The leading terms are removed afterwards.

  begin
    for i in p'range loop
      Clear(ta(i));
      Pop_First_Term(p(i),ta(i));
    end loop;
  end Leading_Terms;

  procedure Find_Max ( ta : in Terms_Array; index : in out Boolean_Array;
                       stop : in out boolean ) is

    res : Degrees := new Natural_Vectors.Vector'(ta(1).dg'range => 0);

  begin
    stop := true;
    for i in ta'range loop
      if ta(i).cf /= CMPLX(0.0)
       then if ta(i).dg > res
             then res.all := Natural_Vectors.Vector(ta(i).dg.all);
                  index(1..(i-1)) := (1..(i-1) => false);
                  index(i) := true;
                  stop := false;
             elsif Equal(ta(i).dg,res)
                 then index(i) := true;
                      stop := false;
            end if;
      end if;
    end loop;
    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(res));
  end Find_Max;

  procedure Update ( p : in out Poly_Sys; n : in natural;
                     ta : in out Terms_Array; da : in out Degrees_Array;
                     nda,cnt : in out natural; mat : in out Matrix;
                     stop : in out boolean ) is

    index : natural;
    is_max : Boolean_Array(1..n) := (1..n => false);

  begin
    Find_Max(ta,is_max,stop);
    if not stop
     then
      -- Get the next Degrees for in the Degrees_Array
       for i in is_max'range loop
         if is_max(i)
          then index := i;  exit;
         end if;
       end loop;
       nda := nda + 1;
       Natural_Vectors.Copy(Natural_Vectors.Link_to_Vector(ta(index).dg),
                            Natural_Vectors.Link_to_Vector(da(nda)));
      -- Fill in the matrix and update the window :
       for i in is_max'range loop
         if is_max(i)
          then mat(i,nda) := ta(i).cf;
               Pop_First_Term(p(i),ta(i));
               cnt := cnt+1;
          else mat(i,nda) := CMPLX(0.0);
         end if;
       end loop;
    end if;
  end Update;

  procedure Coefficient_Matrix
                ( p : in Poly_Sys; mat : in out Matrix;
                  da : in out Degrees_Array; nda : in out natural;
                  diagonal : in out boolean ) is

  -- DESCRIPTION :
  --   Constructs the coefficient matrix of the polynomial system.
  --   Stops when the system is diagonal.

  -- REQUIRED :
  --   mat'range(1) = p'range, mat'range(2) = 1..Sum_Number_of_Terms(p),
  --   da'range = mat'range(2).

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   mat        coefficient matrix, up to column nda filled up;
  --   da         da(1..nda) collects the different terms in the system;
  --   nda        number of different terms;
  --   diagonal   true if the leading terms are all different.

    work : Poly_Sys(p'range);
    stop : boolean := false;
    Window : Terms_Array(p'range);
    cnt : natural := 0;

  begin
    Copy(p,work);
    Leading_Terms(work,Window);
    diagonal := false;
    while not stop and not diagonal loop
      Update(work,p'last,Window,da,nda,cnt,mat,stop);
      if (cnt = p'last) and then (cnt = nda)
       then diagonal := true;
      end if;
      exit when diagonal;
    end loop;
    Clear(work);
  end Coefficient_Matrix;

  procedure Coefficient_Matrix
                ( p : in Poly_Sys; mat : in out Matrix;
                  da : in out Degrees_Array; nda : in out natural ) is

  -- DESCRIPTION :
  --   Constructs the coefficient matrix of the polynomial system.

  -- REQUIRED :
  --   mat'range(1) = p'range, mat'range(2) = 1..Sum_Number_of_Terms(p),
  --   da'range = mat'range(2).

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   mat        coefficient matrix, up to column nda filled up;
  --   da         da(1..nda) collects the different terms in the system;
  --   nda        number of different terms.

    work : Poly_Sys(p'range);
    stop : boolean := false;
    Window : Terms_Array(p'range);
    cnt : natural := 0;

  begin
    Copy(p,work);
    Leading_Terms(work,Window);
    while not stop loop
      Update(work,p'last,Window,da,nda,cnt,mat,stop);
    end loop;
    Clear(work);
  end Coefficient_Matrix;

  procedure Make_Polynomial_System
                    ( p : in out Poly_Sys; mat : in Matrix;
                      da : in Degrees_Array; nda : in natural;
                      inconsistent,infinite : out boolean ) is

    t : Term;
    n : natural := p'length;

  begin
    inconsistent := false;
    infinite := false;
    Clear(p);
    for i in p'range loop
      p(i) := Null_Poly;
      for j in 1..nda loop
        if modulus(mat(i,j)) > mach_eps
         then t.dg := new Natural_Vectors.Vector'(da(j).all); 
              t.cf := mat(i,j);
              Plus_Term(p(i),t);
              Clear(t);
        end if;
      end loop;
      if p(i) = Null_Poly
       then infinite := true;
       elsif Degree(p(i)) = 0
           then inconsistent := true;
      end if; 
    end loop;
  end Make_Polynomial_System;

  function Sum_Number_of_Terms ( p : Poly_Sys ) return natural is

  -- DESCRIPTION :
  --   Returns the sum of the number of terms of every polynomial in p.

    sum : natural := 0;

  begin
    for i in p'range loop
      sum := sum + Number_of_Terms(p(i));
    end loop;
    return sum;
  end Sum_Number_of_Terms;

-- TARGET ROUTINES FOR LINEAR REDUCTION :

  procedure Reduce ( p : in out Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean ) is

    n : natural := p'length;
    cp : Poly_Sys(p'range);
    max_terms : constant natural := Sum_Number_of_Terms(p);
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural := 0;
    mat : Matrix(p'range,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns,diagonal);
    if diagonal
     then inconsistent := false;
          infinite := false;
     else declare
            coeffmat : Matrix(p'range,1..numb_columns);
          begin
            for i in coeffmat'range(1) loop
              for j in coeffmat'range(2) loop
                coeffmat(i,j) := mat(i,j);
              end loop;
            end loop;
            Triangulate(coeffmat,n,numb_Columns);
            Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                                   inconsistent,infinite);
            for i in 1..numb_Columns loop
              Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(Columns(i)));
            end loop;
          end;
    end if;
  end Reduce;

  procedure Sparse_Reduce ( p : in out Poly_Sys;
                            inconsistent,infinite : in out boolean ) is

    n : natural := p'length;
    max_terms : constant natural := Sum_Number_of_Terms(p);
    columns : Degrees_Array(1..max_terms);
    numb_columns : natural := 0;
    mat : Matrix(1..n,1..max_terms);

  begin
    Coefficient_Matrix(p,mat,columns,numb_columns);
    declare
      coeffmat : Matrix(p'range,1..numb_columns);
    begin
      for i in coeffmat'range(1) loop
        for j in coeffmat'range(2) loop
          coeffmat(i,j) := mat(i,j);
        end loop;
      end loop;
      Diagonalize(coeffmat,n,numb_Columns);
      Make_Polynomial_System(p,coeffmat,columns,numb_columns,
                             inconsistent,infinite);
      for i in 1..numb_Columns loop
        Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(columns(i)));
      end loop;
    end;
  end Sparse_Reduce;

-- NONLINEAR REDUCTION :

  function Total_Degree ( p : Poly_Sys ) return natural is

    d : natural := 1;
    tmp : integer;

  begin
    for i in p'range loop
      tmp := Degree(p(i));
      if tmp >= 0
       then d := d * tmp;
      end if;
    end loop;
    return d;
  end Total_Degree;

  function LEQ ( d1,d2 : Degrees ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all degrees of d1 are lower than
  --   or equal to the degrees of d2

  begin
    for i in d1'range loop
      if d1(i) > d2(i)
       then return false;
      end if;
    end loop;
    return true;
  end LEQ;

  function Leading_Term ( p : Poly ) return Term is

  -- DESCRIPTION :
  --   Returns the leading term of the polynomial p.

    tf : Term;

    procedure First_Term (t : in Term; continue : out boolean) is
    begin
      Copy(t,tf);
      continue := false;
    end First_Term;
    procedure Get_First_Term is new Visiting_Iterator (First_Term);
  begin
    Get_First_Term(p);
    return tf;
  end Leading_Term;

  function Can_Be_Eliminated ( p : Poly_Sys; j : natural ) return boolean is

  -- DESCRIPTION :
  --   returns true if the degree of the j-th unknown in each equation 
  --   is zero.

  begin
    for i in p'range loop
      if Degree(p(i),j) > 0
       then return false;
      end if;
    end loop;
    return true;
  end Can_Be_Eliminated;

  procedure Shift_Null_Polynomial ( p : in out Poly_Sys ) is

  -- DESCRIPTION :
  --   The null polynomial in the system p will be shifted down
  --   towards the end.

  begin
    for i in p'range loop
      if p(i) = Null_Poly
       then for j in i..(p'last-1) loop
              Copy(p(j+1),p(j));
              Clear(p(j+1));
            end loop;
      end if;
    end loop;
  end Shift_Null_Polynomial;

  procedure Eliminate ( p : in out Poly; j : in natural ) is

  -- DESCRIPTION :
  --   The j-th unknown will be eliminated out of the polynomial p

    n : natural := Number_Of_Unknowns(p);
    procedure Eliminate_Term (t : in out Term; continue : out boolean) is
      d : Degrees := new Natural_Vectors.Vector(1..(n-1));       
    begin
      for i in 1..(j-1) loop
        d(i) := t.dg(i);
      end loop;
      for i in j..(n-1) loop
        d(i) := t.dg(i+1);
      end loop;
      Clear(t);
      t.dg := d;
      continue := true;
    end Eliminate_Term;
    procedure Eliminate_Terms is new Changing_Iterator(Eliminate_Term);
  begin
    Eliminate_Terms(p);
  end Eliminate;
        
  procedure Eliminate ( p : in out Poly_Sys; j : in natural ) is

  -- DESCRIPTION :
  --   The j-th unknown will be eliminated out of each equation.

  begin
    for i in p'range loop
      Eliminate(p(i),j);
    end loop;
  end Eliminate;
    
  procedure Replace ( p : in out Poly_Sys; pp : in Poly; i : in natural ) is

  -- DESCRIPTION :
  --   This procedure replaces the i-th polynomial in the system p
  --   by the polynomial pp.  If pp is a null polynomial then the procedure
  --   tries to eliminate an unknown, in order to have as much equations
  --   as there are unknowns.

    tmp : natural;

  begin
    if (pp = Null_Poly) or else (Number_Of_Unknowns(pp) = 0)
     then -- try to eliminate an unknown
          tmp := Number_Of_Unknowns(p(1));
          Clear(p(i)); p(i) := Null_Poly;
          for j in reverse 1..Number_Of_Unknowns(p(1)) loop
            if Can_Be_Eliminated(p,j)
             then Eliminate(p,j);
            end if;
          end loop;
          Shift_Null_Polynomial(p);
     else Clear(p(i)); Copy(pp,p(i));
    end if;
  end Replace;

  function red ( p,b1,b2 : Poly ) return Poly is

    Rpb1 : Poly := Rpoly(p,b1);

  begin
    if Number_Of_Unknowns(Rpb1) = 0
     then return Null_Poly;
     else declare
            Rpb2 : Poly := Rpoly(Rpb1,b2);
          begin
            Clear(Rpb1);
            return Rpb2;
          end;
    end if;
  end red;

  function Reduce ( p,b1,b2 : Poly ) return Poly is

  -- DESCRIPTION :
  --   returns p mod < b1,b2 >

    temp : Poly := red(p,b1,b2);
  begin
    if Number_Of_Unknowns(temp) = 0
     then return Null_Poly;
     else Clear(temp);
          return red(p,b2,b1);
    end if;
  end Reduce;

  function Simple_Criterium ( p1,p2 : Poly ) return boolean is

  -- DESCRIPTION :
  --   returns true if lcm(in(p1),in(p2)) = in(p1), if in(p2) | in(p1).

    lt1,lt2 : Term;
    res : boolean;

  begin
    lt1 := Leading_Term(p1);
    lt2 := Leading_Term(p2);
    res := LEQ(lt2.dg,lt1.dg);
    Clear(lt1); Clear(lt2);
    return res;
  end Simple_Criterium;

  procedure Rpoly_Criterium ( p,b1,b2 : in Poly; cnt : in out natural;
                              res : out boolean ) is

  -- DESCRIPTION :
  --   Applies the R-polynomial criterium and counts the number of 
  --   R-polynomials computed.

    Rpb1 : Poly := Rpoly(p,b1);
    Rpb2 : Poly;

  begin
    cnt := cnt + 1;
    if Number_Of_Unknowns(Rpb1) = 0
     then res := true;
     else Rpb2 := Rpoly(Rpb1,b2);
          cnt := cnt + 1;
          Clear(Rpb1);
          if Number_of_Unknowns(Rpb2) = 0
           then res := true;
           else Clear(Rpb2);
                Rpb2 := Rpoly(p,b2);
                cnt := cnt + 1;
                if Number_of_Unknowns(Rpb2) = 0
                 then res := true;
                 else Rpb1 := Rpoly(Rpb2,b1);
                      cnt := cnt + 1;
                      Clear(Rpb2);
                      if Number_of_Unknowns(Rpb1) = 0
                       then res := true;
                       else res := false;
                            Clear(Rpb1);
                      end if;
                end if;
          end if;
    end if;
  end Rpoly_Criterium;

  function Criterium ( p,q,s : Poly ) return boolean is

  -- DESCRIPTION :
  --   returns true if p may be replaced by s.

  begin
    if Simple_Criterium(p,q)
     then return true;
     else declare
            temp : Poly := Reduce(p,q,s);
            res : boolean := (Number_Of_Unknowns(temp) = 0);
          begin
            Clear(temp);
            return res;
          end;
    end if;
  end Criterium;

  procedure Criterium ( p,q,s : in Poly; cnt : in out natural;
                        res : out boolean ) is

  -- DESCRIPTION :
  --   returns true if p may be replaced by s.

  begin
    if Simple_Criterium(p,q)
     then res := true;
     else Rpoly_Criterium(p,q,s,cnt,res);
    end if;
  end Criterium;

  procedure Reduce ( p : in Poly_Sys; res : in out Poly_Sys;
                     cnt_eq : in out natural; max_eq : in natural;
                     cnt_sp : in out natural; max_sp : in natural;
                     cnt_rp : in out natural; max_rp : in natural ) is

    S : Poly;
    n : natural := p'last - p'first + 1;
    dS,dpi,dpj : integer;
    ok : boolean;

    procedure try ( i,dpi : in natural ) is

    -- DESCRIPTION : try to replace p_i by S

      p_red : Poly_Sys(1..n);

    begin
      if cnt_eq > max_eq then return; end if;
      if cnt_sp > max_sp then return; end if;
      Clear(p_red); Copy(p,p_red);
      Replace(p_red,S,i);
      if dS = 0
       then return;
       elsif Total_Degree(p_red) < Total_Degree(res)
           then Copy(p_red,res); 
                Reduce(p_red,res,cnt_eq,max_eq,cnt_sp,max_sp,cnt_rp,max_rp);
           elsif cnt_eq <= max_eq
               then cnt_eq := cnt_eq + 1;
                    Reduce(p_red,res,
                           cnt_eq,max_eq,cnt_sp,max_sp,cnt_rp,max_rp);
      end if;
      Clear(p_red);
    end try;

  begin
    if cnt_eq > max_eq then return; end if;
    if cnt_sp > max_sp then return; end if;
    if cnt_rp > max_rp then return; end if;
    for i in 1..n loop
      for j in (i+1)..n loop
        if (p(i) /= Null_Poly) and (p(j) /= Null_Poly)
         then
           Clear(S); S := Spoly(p(i),p(j));
           cnt_sp := cnt_sp + 1;
           dS  := Degree(S); dpi := Degree(p(i)); dpj := Degree(p(j));
           if dS <= dpi and then dpi > dpj
             and then Criterium(p(i),p(j),S)
            then try(i,dpi);
            elsif dS <= dpj and then dpi < dpj
                 and then Criterium(p(j),p(i),S)
                then try(j,dpj);
                else -- dpi = dpj
                  if dS <= dpi
                   then Criterium(p(i),p(j),S,cnt_rp,ok);
                        if ok then try(i,dpi); end if;
                  end if;
                  if dS <= dpj
                   then Criterium(p(j),p(i),S,cnt_rp,ok);
                        if ok then try(j,dpj); end if;
                  end if;
           end if;
           Clear(S);
        end if;
        exit when (dS = 0);
      end loop;
    end loop;
  end Reduce;

  procedure Sparse_Reduce ( p : in Poly_Sys; res : in out Poly_Sys;
                            cnt_eq : in out natural; max_eq : in natural ) is

    S : Poly;
    n : natural := p'last - p'first + 1;
    dS,dpi,dpj : integer;

    procedure try ( i,dpi : in natural ) is

    -- DESCRIPTION : try to replace p_i by S

      p_red : Poly_Sys(1..n);
      inconsistent,infinite : boolean := false;
    begin
      if cnt_eq > max_eq then return; end if;
      Clear(p_red); Copy(p,p_red);
      Replace(p_red,S,i);
      if dS /= 0
       then Sparse_Reduce(p_red,inconsistent,infinite);
      end if;
      if dS = 0 or inconsistent
       then cnt_eq := max_eq + 1;
            return;
       elsif Total_Degree(p_red) < Total_Degree(res)
           then Copy(p_red,res); 
                Sparse_Reduce(p_red,res,cnt_eq,max_eq);
           else cnt_eq := cnt_eq + 1;
                Sparse_Reduce(p_red,res,cnt_eq,max_eq);
      end if;
      Clear(p_red);
    end try;

  begin
    if cnt_eq > max_eq then return; end if;
    for i in 1..n loop
      for j in (i+1)..n loop
        if (p(i) /= Null_Poly) and (p(j) /= Null_Poly)
         then

              Clear(S); S := Spoly(p(i),p(j));

              dS  := Degree(S);
              dpi := Degree(p(i));
              dpj := Degree(p(j));

              if dS <= dpi and then dpi > dpj
                and then Criterium(p(i),p(j),S)
               then try(i,dpi);
               elsif dS <= dpj and then dpi < dpj
                    and then Criterium(p(j),p(i),S)
                   then try(j,dpj);
                  else -- dpi = dpj
                      if dS <= dpi
                       and then Criterium(p(i),p(j),S) then try(i,dpi); end if;
                      if dS <= dpj
                       and then Criterium(p(j),p(i),S) then try(j,dpj); end if;
              end if;

              Clear(S);

        end if;
      end loop;
    end loop;
  end Sparse_Reduce;

end Reduction_of_Polynomial_Systems;
