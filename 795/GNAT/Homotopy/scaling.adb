with Mathematical_Functions;           use Mathematical_Functions;
with Complex_Numbers;                  use Complex_Numbers;
with Natural_Vectors,Integer_Vectors; 
with Complex_Matrices;                 use Complex_Matrices;
with Complex_Linear_System_Solvers;    use Complex_Linear_System_Solvers;

package body Scaling is

  function log ( b : in natural; x : in double_float ) return double_float is
  begin
    return ( LOG10(x)/LOG10(double_float(b)) );
  end log;

  procedure Scale ( p : in out Poly ) is

    sum : double_float := 0.0;
    number : natural := 0;
    factor : double_complex;

    procedure Add_To_Sum ( t : in Term; continue : out boolean ) is
    begin
      sum := sum + modulus(t.cf);
      number := number + 1;
      continue := true;
    end Add_To_Sum;
    procedure Compute_Sum is new Visiting_Iterator(Add_To_Sum);

  begin
    Compute_Sum(p);
    factor := CMPLX(double_float(number)/sum);
    Mult_Coeff(p,factor);
  end Scale;

  procedure Scale ( s : in out Poly_Sys ) is
  begin
    for i in s'range loop
      scale(s(i));
    end loop;
  end Scale;

 procedure Scale ( s : in out Poly_Sys; bas : in natural := 2;
                   diff : in boolean; cond : out double_float;
                   sccff : out vector ) is
  
  r1,r2,target : Poly;
  n : natural := s'last - s'first + 1;
  nm : natural := 2 * n;
  mat : matrix(1..nm,1..nm);
  right,scaleringscoeff : vector(1..nm);
 
  function center_coefficients ( s : in Poly_Sys; bas : in natural )
                               return Poly is
    r1,r1i : Poly;
    n : natural := s'last - s'first + 1;

    procedure Scan ( p : in Poly; i : in natural ) is

      init : Poly;
      t_init : Term;
   
      procedure Scan_Term ( t : in Term; continue : out boolean ) is
        tt : Term;
        temp : Poly;
      begin
        Copy(init,temp);
        continue := true;
        for k in t.dg'range loop
          if t.dg(k) /= 0
           then tt.cf := CMPLX(double_float(t.dg(k)));
                tt.dg := new Natural_Vectors.Vector'(1..2*n => 0);
                tt.dg(k) := 1;
                Plus_Term(temp,tt);
                Clear(tt);
          end if;
        end loop;
        tt.cf := CMPLX(log(bas,modulus(t.cf)));
        tt.dg := new Natural_Vectors.Vector'(1..2*n => 0);
        Plus_Term(temp,tt);
        Clear(tt);
        Mult_Poly(temp,temp);
        Plus_Poly(r1i,temp);
        Clear(temp);
      end Scan_Term;
      procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

    begin
      t_init.cf := CMPLX(1.0);
      t_init.dg := new Natural_Vectors.Vector'(1..2*n => 0);
      t_init.dg(n+i) := 1;
      init := Create(t_init);
      Clear(t_init);
      Scan_Terms(p);
    end Scan;

  begin
    r1 := Null_Poly;
    for i in s'range loop
      Scan(s(i),i);
      Plus_Poly(r1,r1i);
      Clear(r1i);
    end loop;
    return r1;
  end center_coefficients;

  function reduce_diff ( s : in Poly_Sys; bas : in natural ) return Poly is

    r2,r2i : Poly;

    procedure Scan2 (p : in Poly; t : in Term; nr : in natural) is
      count : natural := 0;
      procedure Scan2_Term (t2 : in Term; continue : out boolean) is
        tt : Term;
        temp : Poly := Null_Poly;
      begin
        continue := true;
        count := count + 1;
        if count > nr
         then for i in t2.dg'range loop
                if t.dg(i) /= t2.dg(i)
                 then tt.dg := new Natural_Vectors.Vector'(1..2*n => 0);
                      tt.dg(i) := 1;
                      tt.cf := CMPLX(double_float(t.dg(i)-t2.dg(i)));
                      Plus_Term(temp,tt);
                      Clear(tt);
                end if;
              end loop;
        end if;
        tt.dg := new Natural_Vectors.Vector'(1..2*n => 0);
        tt.cf := CMPLX(log(bas,(modulus(t.cf)/modulus(t2.cf))));
        Plus_Term(temp,tt);
        Clear(tt);
        Mult_Poly(temp,temp);
        Plus_Poly(r2i,temp);
        Clear(temp);
      end Scan2_Term;
      procedure Scan2_Terms is new Visiting_Iterator(Scan2_Term);
    begin
      Scan2_Terms(p);
    end Scan2;

    procedure Scan ( p : in Poly ) is
      nr : natural := 0;
      procedure Scan_Term ( t : in Term; continue : out boolean ) is
      begin
        nr := nr + 1;
        continue := true;
        Scan2(p,t,nr);
      end Scan_Term;
      procedure Scan_Terms is new Visiting_Iterator(Scan_Term);
    begin
      Scan_Terms(p);
    end Scan;

  begin
    r2 := Null_Poly;
    for i in s'range loop
      Scan(s(i));
      Plus_Poly(r1,r2i);
      Clear(r2i);
    end loop;
    return r2;
  end reduce_diff;

  procedure Make_Linear_System ( r : in Poly; mat : out matrix;
                                 right : out vector ) is
    drj : Poly;

    procedure Init_Linear_System (m : out matrix; r : out vector) is
    begin
      for i in m'range(1) loop
        r(i) := CMPLX(0.0);
        for j in m'range(2) loop
          m(i,j) := CMPLX(0.0);
        end loop;
      end loop;
    end Init_Linear_System;

    procedure Scan (p : in Poly; j : in natural) is

      procedure Scan_Term (t : in Term; continue : out boolean) is
      begin
        continue := true;
        for i in t.dg'range loop
          if t.dg(i) = 1
           then mat(j,i) := t.cf;
                return;
          end if;
        end loop;
        right(j) := -t.cf;
      end Scan_Term;
      procedure Scan_Terms is new Visiting_Iterator (Scan_Term);

    begin
      Scan_Terms(p);
    end Scan;

  begin
   Init_Linear_System(mat,right);
   for j in 1..Number_Of_Unknowns(r) loop
     drj := Complex_Multivariate_Polynomials.Diff(r,j);
     Scan(drj,j);
   end loop;
  end Make_Linear_System;

  procedure Scale ( s : in out Poly_Sys; bas : in natural;
                    mat : in out matrix; right : in out vector;
                    cond : out double_float ) is
 
   n : natural := s'last - s'first + 1;
   ipvt : Integer_Vectors.Vector(1..2*n);
 
   procedure Scale ( p : in out Poly; ip : in natural;
                     scalingcoeff : in vector; bas : in natural ) is

     procedure Scale_Term ( t : in out Term; continue : out boolean ) is
       exp : double_float := 0.0;
     begin
       exp := REAL_PART(scalingcoeff(n+ip));
       for k in t.dg'range loop
         exp := exp + double_float(t.dg(k))*REAL_PART(scalingcoeff(k));
       end loop;
       t.cf := t.cf * CMPLX(double_float(bas) ** exp);
       continue := true;
     end Scale_Term;
     procedure Scale_Terms is new Changing_Iterator (Scale_Term);

   begin
     Scale_Terms(p);
   end scale;
 
  begin
    lufco(mat,2*n,ipvt,cond);
    lusolve(mat,2*n,ipvt,right);
    for i in s'range loop
      scale(s(i),i,right,bas);
    end loop;
  end scale;
        
  begin
    r1 := center_coefficients(s,bas);
    if diff
     then r2 := reduce_diff(s,bas);
          target := r1 + r2;
          clear(r1); clear(r2);
     else copy(r1,target); clear(r1);
    end if;
    Make_Linear_System(target,mat,right);
    clear(target);
    scale(s,bas,mat,right,cond);
    sccff := right;
  end Scale;

  procedure Scale ( basis : in natural; sccff : in Vector;
                    s : in out Solution ) is
  begin
    for i in s.v'range loop
      s.v(i) := CMPLX(double_float(basis)**REAL_PART(sccff(i))) * s.v(i);
    end loop;
  end Scale;

  procedure Scale ( basis : in natural; sccff : in Vector;
                    sols : in out Solution_List ) is
  begin
    if Is_Null(sols)
     then null;
     else declare
            temp : Solution_List := sols;
            n : natural := Head_Of(sols).n;
            s : Solution(n);
            l : Link_To_Solution;
          begin
            while not Is_Null(temp) loop
              l := Head_Of(temp);
              s := l.all;
              Scale(basis,sccff,s);
              Clear(l);
              l := new Solution'(s);
              Set_Head(temp,l);
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
  end Scale;

end Scaling;
