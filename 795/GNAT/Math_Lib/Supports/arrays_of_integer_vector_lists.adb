with unchecked_deallocation;

package body Arrays_of_Integer_Vector_Lists is

-- CONSTRUCTOR :

  procedure Copy ( l1 : in Array_of_Lists; l2 : in out Array_of_Lists ) is
  begin
    for k in l1'range loop
      Copy(l1(k),l2(k));
    end loop;
  end Copy;

-- SELECTORS :

  function Is_Equal ( l1,l2 : Array_of_Lists ) return boolean is
  begin
    if l1'first /= l2'first or else l1'last /= l2'last
     then return false;
     else for k in l1'range loop
            if not Is_Equal(l1(k),l2(k))
             then return false;
            end if;
          end loop;
          return true;
    end if;
  end Is_Equal;

  function Is_Equal ( l1,l2 : Link_to_Array_of_Lists ) return boolean is
  begin
    if l1 = null and then l2 /= null
     then return false;
     elsif l2 = null
         then return true;
         else return Is_Equal(l1.all,l2.all);
    end if;
  end Is_Equal;

  function Length_Of ( l : Array_of_Lists ) return natural is

    res : natural := 0;
 
  begin
    for i in l'range loop
      res := res + Length_Of(l(i));
    end loop;
    return res;
  end Length_Of;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation
       (Array_of_Lists,Link_to_Array_of_Lists);

  procedure Deep_Clear ( l : in out Array_of_Lists ) is
  begin
    for k in l'range loop
      Deep_Clear(l(k));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( l : in out Array_of_Lists ) is
  begin
    for k in l'range loop
      Shallow_Clear(l(k));
    end loop;
  end Shallow_Clear;

  procedure Deep_Clear ( l : in out Link_to_Array_of_Lists ) is
  begin
    if l /= null
     then Deep_Clear(l.all); free(l);
    end if;
  end Deep_Clear;

  procedure Shallow_Clear ( l : in out Link_to_Array_of_Lists ) is
  begin
    if l /= null
     then Shallow_Clear(l.all); free(l);
    end if;
  end Shallow_Clear;

end Arrays_of_Integer_Vector_Lists;
