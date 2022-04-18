with integer_io,Integer_Vectors_io;      use integer_io,Integer_Vectors_io;
with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Complex_Numbers;                    use Complex_Numbers;
with Complex_Multivariate_Laurent_Polynomials;

with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Laurent_to_Polynomial_Converters;   use Laurent_to_Polynomial_Converters;
with Laurent_Polynomial_Randomizers;     use Laurent_Polynomial_Randomizers;

with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Power_Lists;                        use Power_Lists;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
 
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Fewnomials;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;

with Symmetric_BKK_Bound_Solvers;        use Symmetric_BKK_Bound_Solvers;
with Orbits_of_Solutions;                use Orbits_of_Solutions;

package body Symmetric_Polyhedral_Continuation is

  function Symmetric_Mixed_Solve
                ( file : file_type; grp : List_of_Permutations; sign : boolean;
                  p : Laur_Sys; mixsub : Mixed_Subdivision;
                  n : natural; mix : Vector ) return Solution_List is

    sols,sols_last : Solution_List;
    cnt : natural := 0;
    tmp : Mixed_Subdivision := mixsub;

    procedure Solve_Subsystem ( mic : in Mixed_Cell ) is
  
      q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
      sq : Laur_Sys(q'range);
      qsols : Solution_List;
      fail : boolean;
      eps : constant double_float := 10.0**(-10);

    begin
      new_line(file);
      put(file,"*** CONSIDERING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Reduce(n+1,q); sq := Shift(q);
      declare
        pq : Poly_Sys(q'range) := Laurent_to_Polynomial_System(sq);
      begin
        Fewnomials.Solve(sq,qsols,fail);
        if not fail
         then put_line(file,"It is a fewnomial system.");
         else put_line(file,"No fewnomial system.");
              if mic.sub = null
               then put_line(file,"Calling the black box solver.");
                    qsols := Symmetric_BKK_Solve(file,pq,grp,sign);
               else put_line(file,"Using the refinement of the cell.");
                    declare
                      sup : Array_of_Lists(q'range);
                      cnt : natural := sup'first;
                      lif : Array_of_Lists(mix'range);
                      lifq : Laur_Sys(q'range);
                    begin
                      for i in mic.pts'range loop
                        sup(cnt) := Reduce(mic.pts(i),q'last+1);
                        for j in 1..(mix(i)-1) loop
                          Copy(sup(cnt),sup(cnt+j));
                        end loop;
                        cnt := cnt + mix(i);
                      end loop;
                      lif := Induced_Lifting(n,mix,sup,mic.sub.all);
                      lifq := Perform_Lifting(n,mix,lif,q);
                      qsols := Symmetric_Mixed_Solve
                                 (file,grp,sign,lifq,mic.sub.all,n,mix);
                      Deep_Clear(sup); Deep_Clear(lif); Clear(lifq);
                    end;
              end if;
              Set_Continuation_Parameter(qsols,CMPLX(0.0));
        end if;
        put(file,Length_Of(qsols),1);
        put_line(file," solutions found.");
        if not Is_Null(qsols)
         then Analyze(grp,sign,eps,qsols);
              put(file,Length_Of(qsols),1);
              put_line(file," generating solutions found.");
              Mixed_Continuation(file,p,mic.nor.all,qsols);
              Concat(sols,sols_last,qsols);
        end if;
        Clear(pq); Clear(sq);
      end;
      Clear(q); -- Shallow_Clear(qsols);
    end Solve_Subsystem;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      Solve_Subsystem(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return sols;
  end Symmetric_Mixed_Solve;

  function Symmetric_Mixed_Solve
                ( file : file_type; sign : boolean; p : Laur_Sys;
                  mixsub : Mixed_Subdivision; n : natural;
                  mix : Vector ) return Solution_List is

    sols,sols_last : Solution_List;
    cnt : natural;
    tmp : Mixed_Subdivision := mixsub;

    procedure Solve_Subsystem ( mic : in Mixed_Cell ) is
  
      q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
      sq : Laur_Sys(q'range);
      qsols,genqsols : Solution_List;
      fail : boolean;
      eps : constant double_float := 10.0**(-10);

    begin
      new_line(file);
      put(file,"*** CONSIDERING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Reduce(n+1,q); sq := Shift(q);
      declare
        pq : Poly_Sys(q'range) := Laurent_to_Polynomial_System(sq);
      begin
        Fewnomials.Solve(sq,qsols,fail);
        if not fail
         then put_line(file,"It is a fewnomial system.");
         else put_line(file,"No fewnomial system.");
              if mic.sub = null
               then put_line(file,"Calling the black box solver.");
                    qsols := Symmetric_BKK_Solve(file,pq,sign);
               else put_line(file,"Using the refinement of the cell.");
                    declare
                      sup : Array_of_Lists(q'range);
                      cnt : natural := sup'first;
                      lif : Array_of_Lists(mix'range);
                      lifq : Laur_Sys(q'range);
                    begin
                      for i in mic.pts'range loop
                        sup(cnt) := Reduce(mic.pts(i),q'last+1);
                        for j in 1..(mix(i)-1) loop
                          Copy(sup(cnt),sup(cnt+j));
                        end loop;
                        cnt := cnt + mix(i);
                      end loop;
                      lif := Induced_Lifting(n,mix,sup,mic.sub.all);
                      lifq := Perform_Lifting(n,mix,lif,q);
                      qsols := Symmetric_Mixed_Solve(file,sign,lifq,
                                                     mic.sub.all,n,mix);
                      Deep_Clear(sup); Deep_Clear(lif); Clear(lifq);
                    end;
              end if;
              Set_Continuation_Parameter(qsols,CMPLX(0.0));
        end if;
        put(file,Length_Of(qsols),1);
        put_line(file," solutions found.");
        if not Is_Null(qsols)
         then genqsols := Generating(qsols,sign,eps);
              put(file,Length_Of(genqsols),1);
              put_line(file," generating solutions found.");
              Mixed_Continuation(file,p,mic.nor.all,genqsols);
              Concat(sols,sols_last,genqsols);
        end if;
        Clear(pq); Clear(sq);
      end;
      Clear(q); -- Shallow_Clear(genqsols);
    end Solve_Subsystem;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      Solve_Subsystem(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return sols;
  end Symmetric_Mixed_Solve;

end Symmetric_Polyhedral_Continuation;
