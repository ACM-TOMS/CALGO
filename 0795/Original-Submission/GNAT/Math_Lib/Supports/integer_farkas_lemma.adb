--with text_io;                        use text_io;
--with Float_Vectors_io;               use Float_Vectors_io;
--with Float_Matrices_io;              use Float_Matrices_io;

with Floating_Point_Numbers;         use Floating_Point_Numbers;
with Float_Vectors,Float_Matrices;
with Farkas_Lemma;                   use Farkas_Lemma;

package body Integer_Farkas_Lemma is

  procedure Integer_Complementary_Slackness
                  ( tableau : in out matrix; feasible : out boolean ) is
  begin
    Integer_Complementary_Slackness(tableau,tableau'last(2)-1,feasible);
  end Integer_Complementary_Slackness;

  procedure Integer_Complementary_Slackness
                  ( tableau : in out matrix; lastcol : in integer;
                    feasible : out boolean ) is

    tab : Float_Matrices.Matrix
                 (tableau'range(1),tableau'first(2)..lastcol);
    rhs,sol : Float_Vectors.Vector(tab'range(1));
   -- tol : constant single_float := 10.0**(-5);  -- single precision
    tol : constant double_float := 10.0**(-12);  -- double precision
    columns : vector(sol'range);

  begin
    for i in tab'range(1) loop
      for j in tab'range(2) loop
        tab(i,j) := double_float(tableau(i,j));
      end loop;
    end loop;
    for i in rhs'range loop
      rhs(i) := double_float(tableau(i,tableau'last(2)));
    end loop;
   -- put_line("The tableau : "); put(tab,3,3,3);
   -- put_line(" with right hand side : "); put(rhs,3,3,3); new_line;
    Complementary_Slackness(tab,lastcol,rhs,tol,sol,columns,feasible);
   -- put_line("The solution : "); put(sol,3,3,3); new_line;
  end Integer_Complementary_Slackness;

end Integer_Farkas_Lemma;
