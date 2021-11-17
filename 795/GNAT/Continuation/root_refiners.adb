with integer_io,Complex_Numbers_io;      use integer_io,Complex_Numbers_io;
with Complex_Numbers,Complex_Norms;      use Complex_Numbers,Complex_Norms;
with Complex_Matrices,Integer_Vectors;   use Complex_Matrices;
with Complex_Linear_System_Solvers;      use Complex_Linear_System_Solvers;
with Floating_Equalities;                use Floating_Equalities;
with Solutions_io;                       use Solutions_io;

package body Root_Refiners is

  use Floating_Point_Numbers.double_float_io;

-- AUXILIARIES :

  function Is_Real ( sol : Solution; tol : double_float ) return boolean is
  begin
    for i in sol.v'range loop
      if ABS(IMAG_PART(sol.v(i))) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Real;

  function Is_Equal ( s1,s2 : Solution; tol : double_float ) return boolean is
  begin
    for i in s1.v'range loop
      if not Is_Equal(s1.v(i),s2.v(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Complex_Conjugate ( s : Solution ) return Solution is

    res : Solution(s.n) := s;

  begin
    for i in res.v'range loop
      res.v(i) := CMPLX(REAL_PART(s.v(i)),-IMAG_PART(s.v(i)));
    end loop;
    return res;
  end Complex_Conjugate;

  function Is_Clustered ( sol : Solution; nb : natural; sols : Solution_List;
                          tol : double_float ) return natural is

    tmp : Solution_List := sols;
    cnt : natural := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt /= nb
       then if Is_Equal(sol,Head_Of(tmp).all,tol)
	     then return cnt;
            end if;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return nb;
  end Is_Clustered;

  function Is_Clustered ( sol : Solution; nb : natural; sols : Solution_Array;
                          tol : double_float ) return natural is
  begin
    for i in sols'range loop
      if i /= nb
       then if Is_Equal(sol,sols(i).all,tol)
             then return i;
            end if;
      end if;
    end loop;
    return nb;
  end Is_Clustered;

  function Multiplicity ( sol : Solution; sols : Solution_List; 
                          tol : double_float ) return natural is

    tmp : Solution_List := sols;
    cnt : natural := 0;

  begin
    while not Is_Null(tmp) loop
      if Is_Equal(sol,Head_Of(tmp).all,tol)
       then cnt := cnt + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return cnt;
  end Multiplicity;

  function Multiplicity ( sol : Solution; sols : Solution_Array;
                          tol : double_float ) return natural is
    cnt : natural := 0;

  begin
    for i in sols'range loop
      if Is_Equal(sol,sols(i).all,tol)
       then cnt := cnt + 1;
      end if;
    end loop;
    return cnt;
  end Multiplicity;

  procedure Write_Bar ( file : in file_type ) is
  begin
    for i in 1..75 loop
      put(file,'=');
    end loop;
    new_line(file);
  end Write_Bar;

  procedure Write_Info ( file : in file_type; zero : in Solution;
                         i,n,numb : in natural; fail : in boolean ) is

  -- DESCRIPTION :
  --   The information concerning the zero is written

  begin
   -- Write_Bar(file);
    put(file,"solution : "); put(file,i,1); put(file," :        ");
    put(file," start residual : "); put(file,zero.res,2,3,3); new_line(file);
    put(file,zero);
  end Write_Info;
 
  procedure Root_Accounting
               ( file : in file_type; ls : in out Link_to_Solution;
                 nb : in natural; sa : in out Solution_Array;
                 fail : in boolean; tolsing,tolclus : in double_float;
                 nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus : in out natural ) is

  -- DESCRIPTION :
  --   This procedure does root accounting of the solution sol, w.r.t. the
  --   solution list sols.  Information will be provided concerning the type
  --   of solution.

  begin
    if fail
     then put_line(file," no solution ==");
          nbfail := nbfail + 1;
          ls.m := 0;
     else if Is_Real(ls.all,10.0**(-13))
	   then put(file," real ");
	        nbreal := nbreal + 1;
	   else put(file," complex ");
	        nbcomp := nbcomp + 1;
          end if;
          if sa(nb).rco < tolsing
	   then declare
                  m : natural := Multiplicity(ls.all,sa,tolclus);
                begin
                  if m = 1
                   then m := 0;
		  end if;
		  ls.m := m;
	        end;
	        put_line(file,"singular ==");
	        nbsing := nbsing + 1;
           else declare
                  nb2 : natural := Is_Clustered(ls.all,nb,sa,tolclus);
                begin
                  if nb2 = nb
		   then put_line(file,"regular ==");
		        nbreg := nbreg + 1;
                   else put(file,"clustered : ");
		        put(file,nb2,1);
		        put_line(file," ==");
		        nbclus := nbclus + 1;
                  end if;
		  ls.m := 1;
	        end;	   
          end if;  
    end if;
  end Root_Accounting;

  procedure Root_Accounting 
                ( ls : in out Link_to_Solution; nb : in natural;
                  sa : in out Solution_Array; fail : in boolean;
                  tolsing,tolclus : in double_float ) is

  -- DESCRIPTION :
  --   This procedure does root accounting of the solution sol, w.r.t. the
  --   solution list sols.  Information will be provided concerning the type
  --   of solution.

  begin
    if fail
     then ls.m := 0;
     elsif sa(nb).rco < tolsing
         then declare
                m : natural := Multiplicity(ls.all,sa,tolclus);
              begin
                if m = 1
                 then ls.m := 0;
                 else ls.m := m;
                end if;
              end;
         else ls.m := 1;
    end if;
  end Root_Accounting;

  procedure Write_Global_Info
             ( file : in file_type;
               tot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus : in natural ) is

  begin
    Write_Bar(file);
    put(file,"A list of "); put(file,tot,1);
    put_line(file," solutions has been refined :");
    put(file,"Number of regular solutions   : "); put(file,nbreg,1);
    put_line(file,".");
    put(file,"Number of singular solutions  : "); put(file,nbsing,1);
    put_line(file,".");
    put(file,"Number of real solutions      : "); put(file,nbreal,1);
    put_line(file,".");
    put(file,"Number of complex solutions   : "); put(file,nbcomp,1);
    put_line(file,".");
    put(file,"Number of clustered solutions : "); put(file,nbclus,1);
    put_line(file,".");
    put(file,"Number of failures            : "); put(file,nbfail,1);
    put_line(file,".");
    Write_Bar(file);
  end Write_Global_Info;

-- TARGET ROUTINES :

  procedure Silent_Newton
               ( p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jacobi;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural; max : in natural;
                 fail : out boolean ) is

    n : natural := p_eval'length;
    jacobian : matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    y,deltax : vector(1..n);

  begin
    y := eval(p_eval,zero.v);               -- y = f(zero)
    for i in 1..max loop
      jacobian := eval(j_eval,zero.v);      -- solve jacobian*deltax = -f(zero)
     -- Scale(jacobian,y);
      lufco(jacobian,n,ipvt,zero.rco);
      if 1.0 + zero.rco = 1.0
       then fail := (norm1(y) > epsfa);
            return;                         -- singular Jacobi matrix
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);  
      Plus_Vector(zero.v,deltax);           -- make the updates
      y := eval(p_eval,zero.v);
      zero.err := norm1(deltax); zero.res := norm1(y);
      numit := numit + 1;
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) 
                                            -- stopping criteria
       then fail := false; exit;
       elsif numit >= max
           then fail := true; exit;
      end if;
    end loop;
    jacobian := eval(j_eval,zero.v);        -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
    when numeric_error | constraint_error => fail := true; return;
  end Silent_Newton;

  procedure Reporting_Newton
               ( file : in file_type;
                 p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jacobi;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural; max : in natural;
                 fail : out boolean ) is

    n : natural := p_eval'length;
    jacobian : matrix(1..n,1..n);
    ipvt : Integer_Vectors.Vector(1..n);
    y,deltax : vector(1..n);

  begin
    y := Eval(p_eval,zero.v);              -- y = f(zero)
    for i in 1..max loop
      jacobian := eval(j_eval,zero.v);     -- solve jacobian*deltax = -f(zero)
     -- Scale(jacobian,y);
      lufco(jacobian,n,ipvt,zero.rco);
      if 1.0 + zero.rco = 1.0              -- singular jacobian matrix
       then fail := (norm1(y) > epsfa);    -- accuracy not reached yet
            return;
      end if;
      deltax := -y;
      lusolve(jacobian,n,ipvt,deltax);  
      Plus_Vector(zero.v,deltax);          -- make the updates
      y := eval(p_eval,zero.v);
      zero.err := norm1(deltax); zero.res := norm1(y);
      numit := numit + 1;
      put(file,"Step "); put(file,numit,4); new_line(file);      -- output
      put(file," |errxa| : "); put(file,zero.err); new_line(file);
      put(file," |errfa| : "); put(file,zero.res); new_line(file);
      if ( zero.err < epsxa ) or else ( zero.res < epsfa ) 
                                                  -- stopping criteria
       then fail := false; exit;
       elsif numit >= max
           then fail := true; exit;
      end if;
    end loop;
    jacobian := eval(j_eval,zero.v);              -- compute condition number
    lufco(jacobian,n,ipvt,zero.rco);
  exception
    when numeric_error | constraint_error => fail := true; return;
  end Reporting_Newton;

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural ) is

    n : natural := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jacobi(1..n,1..n) := Create(p);
    jac_eval : Eval_Jacobi(1..n,1..n) := Create(jac);
    numb : natural;
    fail : boolean;
    sa : Solution_Array(1..Length_Of(sols)) := Create(sols);

  begin
    for i in sa'range loop
      numb := 0;
      sa(i).res := Norm1(Eval(p_eval,sa(i).v));
      if sa(i).res < 1.0
       then Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
       else fail := true;
      end if;
      Root_Accounting(sa(i),i,sa(sa'first..i),fail,tolsing,epsxa);
      numit := numit + numb;
    end loop;
    clear(p_eval); clear(jac); clear(jac_eval);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural ) is

    n : natural := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jacobi(1..n,1..n) := Create(p);
    jac_eval : Eval_Jacobi(1..n,1..n) := Create(jac);
    numb : natural;
    fail : boolean;
    sa : Solution_Array(1..Length_Of(sols)) := Create(sols);
    refsols_last : Solution_List;

  begin
    for i in sa'range loop
      numb := 0;
      sa(i).res := Norm1(Eval(p_eval,sa(i).v));
      if sa(i).res < 1.0
       then Silent_Newton(p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
            if not fail
             then Append(refsols,refsols_last,sa(i).all);
            end if;
       else fail := true;
      end if;
      Root_Accounting(sa(i),i,sa(sa'first..i),fail,tolsing,epsxa);
      numit := numit + numb;
    end loop;
    clear(p_eval); clear(jac); clear(jac_eval);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
  end Silent_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural;
                 wout : in boolean ) is

    n : natural := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jacobi(1..n,1..n) := Create(p);
    jac_eval : Eval_Jacobi(1..n,1..n) := Create(jac);
    numb : natural;
    nbfail,nbreg,nbsing,nbclus,nbreal,nbcomp : natural := 0;
    nbtot : constant natural := Length_Of(sols);
    fail : boolean;
    sa : Solution_Array(1..nbtot) := Create(sols);

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    Write_Bar(file);
    for i in sa'range loop 
      numb := 0;
      sa(i).res := Norm1(Eval(p_eval,sa(i).v));
      if sa(i).res < 1.0
       then
         if wout
          then Reporting_Newton
                 (file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
          else Silent_Newton
                 (p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
         end if;
       else 
         fail := true;
      end if;
      Write_Info(file,sa(i).all,i,n,numb,fail);
      Root_Accounting(file,sa(i),i,sa(sa'first..i),fail,tolsing,epsxa,
                      nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
      numit := numit + numb;
    end loop;
    Write_Global_Info(file,nbtot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
    clear(p_eval); clear(jac); clear(jac_eval);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural;
                 wout : in boolean ) is

    n : natural := p'length;
    p_eval : Eval_Poly_Sys(1..n) := Create(p);
    jac : Jacobi(1..n,1..n) := Create(p);
    jac_eval : Eval_Jacobi(1..n,1..n) := Create(jac);
    numb : natural;
    nbfail,nbreg,nbsing,nbclus,nbreal,nbcomp : natural := 0;
    nbtot : constant natural := Length_Of(sols);
    fail : boolean;
    sa : Solution_Array(1..nbtot) := Create(sols);
    refsols_last : Solution_List;

  begin
    new_line(file);
    put_line(file,"THE SOLUTIONS :"); new_line(file);
    put(file,nbtot,1); put(file," "); put(file,n,1); new_line(file);
    Write_Bar(file);
    for i in sa'range loop 
      numb := 0;
      sa(i).res := Norm1(Eval(p_eval,sa(i).v));
      if sa(i).res < 1.0
       then
         if wout
          then Reporting_Newton
                 (file,p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
          else Silent_Newton
                 (p_eval,jac_eval,sa(i).all,epsxa,epsfa,numb,max,fail);
         end if;
         if not fail
          then Append(refsols,refsols_last,sa(i).all);
         end if;
       else 
         fail := true;
      end if;
      Write_Info(file,sa(i).all,i,n,numb,fail);
      Root_Accounting(file,sa(i),i,sa(sa'first..i),fail,tolsing,epsxa,
                      nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
      numit := numit + numb;
    end loop;
    Write_Global_Info(file,nbtot,nbfail,nbreal,nbcomp,nbreg,nbsing,nbclus);
    clear(p_eval); clear(jac); clear(jac_eval);
    Deep_Clear(sols); sols := Create(sa); Clear(sa);
  end Reporting_Root_Refiner;

end Root_Refiners;
