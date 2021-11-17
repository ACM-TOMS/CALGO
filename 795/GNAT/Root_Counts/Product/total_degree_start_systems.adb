with text_io;                           use text_io;
with Random_Number_Generators;          use Random_Number_Generators;
with Complex_Numbers;                   use Complex_Numbers;
with Natural_Vectors;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package body Total_Degree_Start_Systems is

  procedure Total_Degree_Info is

  -- DESCRIPTION :
  --   Displays information about the total degree on screen.

    i : array(1..5) of string(1..65);

  begin
    i(1):="  The  total  degree  is  the  product  of  the  degrees  of  the";
    i(2):="polynomials in the system.  The i-th equation of the start system";
    i(3):="is a univariate polynomial in the i-th unknown of the same degree";
    i(4):="as  the i-th polynomial in the system that has to be solved.  The";
    i(5):="total degree equals the number of solutions of the start system. ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Total_Degree_Info;

  procedure Start_Solutions
               ( level : in natural;
                 q : in Poly_Sys; c : in Vector; s : in out Solution;
                 qsols,qsols_last : in out Solution_List ) is

  -- DESCRIPTION :
  --   All solutions to the polynomial system q are computed.
  --   The parameter level indicates the current component in the recursive
  --   application of the rule of de Moivre.

    d : natural;

  begin
    if level <= s.n
     then d := Degree(q(level));
          for j in 1..d loop
            s.v(level) := de_Moivre(d,j,c(level));
            Start_Solutions(level+1,q,c,s,qsols,qsols_last);
          end loop;
     else s.t := CMPLX(0.0);
          s.m := 1;
          s.err := 0.0; s.rco := 1.0; s.res := 0.0;
          Append(qsols,qsols_last,s);
    end if;
  end Start_Solutions;

  procedure Start_System 
               ( p : in Poly_Sys; q : in out Poly_Sys; c : in Vector;
                 qsols : in out Solution_List ) is
  
    t : Term;
    n : natural := p'length;
    s : Solution(n);
    last : Solution_List := qsols;

  begin
    for i in p'range loop
      t.dg := new Natural_Vectors.Vector'(1..n => 0);
      t.dg(i) := Degree(p(i));
      t.cf := CMPLX(1.0);
      q(i) := Create(t);
      Clear(t);
      t.dg := new Natural_Vectors.Vector'(1..n => 0);
      t.cf := -c(i);
      Plus_Term(q(i),t);
      Clear(t);
    end loop;
    Start_Solutions(1,q,c,s,qsols,last);
  end Start_System;
 
  procedure Start_System
               ( p : in Poly_Sys; q : in out Poly_Sys;
                 qsols : in out Solution_List ) is

    n : natural := p'length;
    c : Vector(1..n);

  begin
    for i in c'range loop
      c(i) := Random1;
    end loop;
    Start_System(p,q,c,qsols);
  end Start_System;

end Total_Degree_Start_Systems;
