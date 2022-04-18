package body Integer_Graded_Lexicographical_Ordening is

  function "<" ( v1,v2 : Vector ) return boolean is

    s1,s2 : integer;
  begin
    s1 := Sum(v1);
    s2 := Sum(v2);
    if s1 < s2
     then return true;
     elsif s1 > s2
         then return false;
         else if v1'first /= v2'first or else v1'last /= v2'last
               then raise Range_Error;
               else for i in v1'range loop
                      if v1(i) < v2(i)
                       then return true;
                       elsif v1(i) > v2(i)
                           then return false;
                      end if;
                    end loop;
                    return false;  -- v1 = v2
              end if;
    end if;
  end "<";

  function "<" ( v1,v2 : Link_to_Vector ) return boolean is
  begin
    if v2 = null
     then return false;
     elsif v1 = null
         then if Sum(v2) > 0
               then return true;
               else return false;
              end if;
         else return v1.all < v2.all;
    end if;
  end "<";

  function ">" ( v1,v2 : Vector ) return boolean is

    s1,s2 : integer;
  begin
    s1 := Sum(v1);
    s2 := Sum(v2);
    if s1 < s2
     then return false;
     elsif s1 > s2
         then return true;
         else if v1'first /= v2'first or else v1'last /= v2'last
               then raise Range_Error;
               else for i in v1'range loop
                      if v1(i) < v2(i)
                       then return false;
                       elsif v1(i) > v2(i)
                           then return true;
                      end if;
                    end loop;
                    return false;  -- v1 = v2
              end if;
    end if;
  end ">";

  function ">" ( v1,v2 : Link_to_Vector ) return boolean is

  begin
    if v1 = null
     then return false;
     elsif v2 = null
         then if Sum(v1) > 0
               then return true;
               else return false;
              end if;
         else return v1.all > v2.all;
    end if;
  end ">";

end Integer_Graded_Lexicographical_Ordening;
