with Greatest_Common_Divisor;         use Greatest_Common_Divisor;
with Integer_Linear_System_Solvers;   use Integer_Linear_System_Solvers;  

--with text_io,integer_io;      use text_io,integer_io;
--with Integer_Vectors_io;      use Integer_Vectors_io;
--with Integer_Matrices_io;     use Integer_Matrices_io;

package body Integer_Linear_Inequality_Solvers is

  procedure Triangulate ( l : in natural; m : in matrix;
                          index : in natural; ineq : in out matrix ) is

    a,b,lcmab,faca,facb : integer;

  begin
    if ineq(index,l) /= 0
     then a := m(l,l); b := ineq(index,l);
          lcmab := lcm(a,b);
          if lcmab < 0 then lcmab := -lcmab; end if;
          faca := lcmab/a;  facb := lcmab/b;
          if facb < 0
           then facb := -facb; faca := -faca;
          end if;
          for j in l..ineq'last(2) loop
            ineq(index,j) := facb*ineq(index,j) - faca*m(l,j);
          end loop;
    end if;
  end Triangulate;

  procedure Triangulate ( l : in natural; m : in matrix;
                          first,last : in natural; ineq : in out matrix ) is
  begin
    for i in first..last loop
      if ineq(i,l) /= 0
       then Triangulate(l,m,i,ineq);
      end if;
    end loop;
  end Triangulate;

  procedure Triangulate ( index,start : in natural; ineq : in out matrix ) is
  
    column : natural := start;  -- current column counter
    firstineq : natural := ineq'first(1);
    found : boolean;
    ind2 : natural;
    a,b,lcmab,faca,facb : integer;

  begin
    --put_line("INEQUALTITY");
    --put(" BEFORE UPDATE : ");
    --for l in ineq'range(2) loop
    --  put(ineq(index,l),1); put(' ');
    --end loop;
    --new_line;
    loop
     -- SEARCH FOR FIRST NONZERO ENTRY IN CURRENT INEQUALITY :
      while column < ineq'last(2) and then ineq(index,column) = 0 loop
        column := column + 1;
      end loop;
      exit when ((ineq(index,column) = 0) or else (column = ineq'last(2)));
                                        -- nothing to eliminate
     -- SEARCH FOR INEQUALITY,
     --    WITH SAME FIRST NONZERO COLUMN, BUT WITH OPPOSITE SIGN :
      found := false;
      for k in firstineq..(index-1) loop
        if ineq(index,column)*ineq(k,column) < 0 -- check for sign
         then found := true;
              for l in start..column-1 loop      -- check if same zero pattern
               -- if ineq(index,l) = 0
               --  then 
                found := (ineq(k,l) = 0);
               -- end if;
                exit when not found;
              end loop;
              if found then ind2 := k; end if;
        end if;
        exit when found;
      end loop;
      exit when not found;  -- no possibility for elimination
     -- if found
     --  then
        -- ELIMINATE BY MAKING A POSITIVE LINEAR COMBINATION :
         a := ineq(index,column);  b := ineq(ind2,column);
         if a < 0
          then lcmab := lcm(-a,b);
               faca := lcmab/(-a); facb := lcmab/b;
          else lcmab := lcm(a,-b);
               facb := lcmab/(-b); faca := lcmab/a;
         end if;
         if ineq(index,ineq'last(2)) >= 0
           or else  -- PRESERVE SIGN OF ineq(index,ineq'last(2)) !!!!
                (faca*ineq(index,ineq'last(2)) 
                            + facb*ineq(ind2,ineq'last(2)) < 0)
          then
            for l in start..ineq'last(2) loop
              ineq(index,l) := faca*ineq(index,l) + facb*ineq(ind2,l);
            end loop;
           -- PROCEED :
            column := column + 1;
            firstineq := ineq'first(1);
          else
           -- TRY TO ELIMINATE WITH OTHER INEQUALITIES :
            firstineq := ind2 + 1;
         end if;
         if (firstineq >= index)
           -- impossible to eliminate with sign preservation
          then firstineq := ineq'first(2);  
               column := column + 1;
         end if;
     --  else
     --    column := column + 1;
     --    firstineq := ineq'first(2);
     -- end if;
      exit when (column >= ineq'last(2)-1);
    end loop;
    --put(" AFTER UPDATE : ");
    --for l in ineq'range(2) loop
    --  put(ineq(index,l),1); put(' ');
    --end loop;
    --new_line;
  end Triangulate;

end Integer_Linear_Inequality_Solvers;
