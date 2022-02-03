with Integer_Support_Functions;          use Integer_Support_Functions;
with Transformations;                    use Transformations;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Lists_of_Vectors_Utilities;         use Lists_of_Vectors_Utilities;

package body Arrays_of_Lists_Utilities is

  function All_Equal ( al : Array_of_Lists ) return boolean is
  begin
    for i in (al'first+1)..al'last loop
      if not Is_Equal(al(al'first),al(i))
       then return false;
      end if;
    end loop;
    return true;
  end All_Equal;

  function Interchange2 ( al : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(al'range);
    index : integer;

  begin
    if Length_Of(al(al'first)) <= 2
     then res := al;
     else index := al'first;
          for i in al'first+1..al'last loop
            if Length_Of(al(i)) <= 2
             then index := i;
             else res(i) := al(i);
            end if;
            exit when index > al'first;
          end loop;
          if index = al'first
           then res(index) := al(index);
           else res(index) := al(al'first);
                res(res'first) := al(index);
                res(index+1..res'last) := al(index+1..al'last);
          end if;
    end if;
    return res;
  end Interchange2;

  function Index2 ( al : Array_of_Lists ) return integer is
  begin
    for i in al'range loop
      if Length_Of(al(i)) <= 2
       then return i;
      end if;
    end loop;
    return al'first;
  end Index2;

  procedure Mixture ( al : in Array_of_Lists;
                      perm,mix : out Link_to_Vector ) is

    wrkper,wrkmix : vector(al'range);    -- intermediate results
    nbd : natural := 0;                  -- # different sets
    ind,min : integer;

    procedure Sort ( indal,indmix : in natural ) is
 
    -- DESCRIPTION :
    --   Puts all lists which are equal to al(perm(index)) together.

    -- ON ENTRY :
    --   indal       the current entry in al;
    --   indmix      the current entry in wrkmix.

    begin
      for j in indal+1..al'last loop
        if Is_Equal(al(wrkper(indal)),al(wrkper(j)))
         then if j /= indal + wrkmix(indmix)
               then declare
                      pos : natural := indal + wrkmix(indmix);
                      tmppos : natural;
                    begin
                      tmppos := wrkper(j);
                      wrkper(j) := wrkper(pos);
                      wrkper(pos) := tmppos;
                    end;
              end if;
              wrkmix(indmix) := wrkmix(indmix) + 1;
        end if;
      end loop;
    end Sort;

    procedure Permute ( ind,nb : in natural ) is

    -- DESCRIPTION :
    --   Changes the permutation vector such that the entry given by
    --   the index stands in front.  The number of different supports is
    --   given by the parameter nb.

      newper : vector(wrkper'range);
      cntnew : natural := newper'first + wrkmix(ind);
      cntwrk : natural := wrkper'first;

    begin
      for i in 1..nb loop
        if i /= ind 
         then for j in 0..wrkmix(i)-1 loop
                newper(cntnew+j) := wrkper(cntwrk+j);
              end loop;
              cntnew := cntnew + wrkmix(i);
         else for j in 0..wrkmix(ind)-1 loop
                newper(newper'first+j) := wrkper(cntwrk+j);
              end loop;
        end if;
        cntwrk := cntwrk + wrkmix(i);
      end loop;
      wrkper := newper;
    end Permute;

  begin
   -- INITIALIZATIONS :
    for i in wrkper'range loop
      wrkper(i) := i;
    end loop;
    wrkmix := (wrkmix'range => 1);
   -- SORTING THE SETS :
    ind := al'first;
    while ind <= al'last loop
      nbd := nbd + 1;
      Sort(ind,nbd);
      ind := ind + wrkmix(nbd);
    end loop;
   -- MINIMAL OCCURENCE SHOULD APPEAR FIRST :
    ind := wrkmix'first;
    min := wrkmix(ind);
    for i in wrkmix'first+1..nbd loop
      if wrkmix(i) < min
       then min := wrkmix(i); ind := i;
      end if;
    end loop;
   -- put("The type of mixture : " ); put(wrkmix(wrkmix'first..nbd)); new_line;
   -- put("The permutation vector : "); put(wrkper); new_line;
    if ind /= wrkmix'first
     then Permute(ind,nbd);
          wrkmix(ind) := wrkmix(wrkmix'first);
          wrkmix(wrkmix'first) := min;
    end if;
   -- put("The type of mixture : " ); put(wrkmix(wrkmix'first..nbd)); new_line;
   -- put("The permutation vector : "); put(wrkper); new_line;
   -- RETURNING THE RESULTS :
    perm := new Integer_Vectors.Vector'(wrkper);
    mix := new Integer_Vectors.Vector'(wrkmix(wrkmix'first..nbd));
  end Mixture;

  function Permute ( perm : Vector; al : in Array_of_Lists )
                   return Array_of_Lists is

    res : Array_of_Lists(al'range);

  begin
    for i in al'range loop
      res(i) := al(perm(i));
    end loop;
    return res;
  end Permute;

  function Different_Points ( al : Array_of_Lists ) return List is

    tmp,res,res_last : List;

  begin
    for i in (al'first+1)..al'last loop
      tmp := al(i);
      while not Is_Null(tmp) loop
        declare
          lv : Link_to_Vector := Head_Of(tmp);
        begin
          if not Is_In(res,lv.all)
           then Append(res,res_last,lv.all);
          end if;
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Different_Points;

  function Different_Points ( al : Array_of_Lists ) return Array_of_Lists is

    res : Array_of_Lists(al'range);

  begin
    res(res'first) := al(al'first);
    for i in (al'first+1)..al'last loop
      res(i) := Different_Points(al(i));
    end loop;
    return res;
  end Different_Points;

  procedure Remove_Duplicates ( al : in out Array_of_Lists ) is
  begin
    for i in al'range loop
      Remove_Duplicates(al(i));
    end loop;
  end Remove_Duplicates;

  procedure Shift ( al : in out Array_of_Lists;
                    shiftvecs : in Integer_Vectors_of_Vectors.Vector ) is
  begin
    for k in al'range loop
      Shift(al(k),shiftvecs(k));
    end loop;
  end Shift;

  function Shift ( al : Array_of_Lists;
                   shiftvecs : Integer_Vectors_of_Vectors.Vector )
                 return Array_of_Lists is

    res : Array_of_Lists(al'range);

  begin
    for k in res'range loop
      res(k) := Shift(al(k),shiftvecs(k));
    end loop;
    return res;
  end Shift;

  procedure Projection ( al : in Array_of_Lists; v : in Vector;
                         ind : integer; res : in out Array_of_Lists;
                         degenerate : out boolean ) is

    pv : integer;
    t : Transfo := Build_Transfo(v,ind);

    procedure Clean ( i : in integer ) is
    begin
      for j in res'first..i loop
        Deep_Clear(res(j));
      end loop;
      Clear(t);
    end Clean;

  begin
    degenerate := false;
    for i in res'range loop
      declare
        pvl : List;
        l : List renames al(i+1);
      begin
        pv := Maximal_Support(l,v);
        pvl := Face(l,v,pv);
        if Length_Of(pvl) <= 1
         then degenerate := true;
              Deep_Clear(pvl);  Clean(i);
              return;
         else res(i) := Transform_and_Reduce(t,ind,pvl);
              Remove_Duplicates(res(i));
              if Length_Of(res(i)) <= 1
               then degenerate := true;
                    Deep_Clear(pvl); Clean(i);
                    return;
              end if;
        end if;
        Deep_Clear(pvl);
      end;
    end loop;
    Clear(t);
  end Projection;

end Arrays_of_Lists_Utilities;
