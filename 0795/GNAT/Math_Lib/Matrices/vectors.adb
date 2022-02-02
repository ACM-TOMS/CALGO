with unchecked_deallocation;

package body Vectors is

  procedure Clear ( v : in out Vector ) is
  begin
    for i in v'range loop
      clear(v(i));
    end loop;
  end Clear;

  function Equal ( v1,v2 : Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;
     else for i in v1'range loop
            if not equal(v1(i),v2(i))
             then return false;
            end if;
          end loop;
          return true;
    end if;
  end Equal;

  procedure Copy ( v1: in Vector; v2 : in out Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else Clear(v2);
          for i in v1'range loop
            copy(v1(i),v2(i));
          end loop;
    end if;
  end Copy;

  function "+" ( v1,v2 : Vector ) return Vector is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else declare
            res : Vector(v1'range);
          begin
            for i in v1'range loop
              res(i) := v1(i) + v2(i);
            end loop;
            return res;
          end;
    end if;
  end "+";
  
  function "-" ( v : Vector ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := -v(i);
    end loop;
    return res;
  end "-";

  function "-" ( v1,v2 : Vector ) return Vector is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else declare
            res : Vector(v1'range);
          begin
            for i in v1'range loop
              res(i) := v1(i) - v2(i);
            end loop;
            return res;
          end;
    end if;
  end "-";

  function "*" ( v : Vector; a : coefftp ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := v(i) * a;
    end loop;
    return res;
  end "*";

  function "*" ( a : coefftp; v : Vector ) return Vector is
  begin
    return v*a;
  end "*";

  function "*" ( v1,v2 : Vector ) return coefftp is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else declare 
            temp,sum : coefftp := zero;
          begin
            for i in v1'range loop
              temp := v1(i)*v2(i);
              Plus_Coeff(sum,temp);
              clear(temp);
            end loop;
            return sum;
          end;
    end if;
  end "*";

  function "*" ( v1,v2 : Vector ) return Vector is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else declare
            res : Vector(v1'range);
          begin
            for i in v1'range loop
              res(i) := v1(i)*v2(i);
            end loop;
            return res;
          end;      
    end if;
  end "*";

  function Sum ( v : Vector ) return coefftp is

    res : coefftp := zero;

  begin
    for i in v'range loop
      Plus_Coeff(res,v(i));
    end loop;
    return res;
  end Sum;

  procedure Plus_Vector ( v1 : in out Vector; v2 : in Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else for i in v1'range loop
            Plus_Coeff(v1(i),v2(i));
          end loop;
    end if;
  end Plus_Vector;

  procedure Min_Vector ( v : in out Vector ) is
  begin
    for i in v'range loop
      Min_Coeff(v(i));
    end loop;
  end Min_Vector;

  procedure Min_Vector ( v1 : in out Vector; v2 : in Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else for i in v1'range loop
            Min_Coeff(v1(i),v2(i));
          end loop;
    end if;
  end Min_Vector;

  procedure Mult_Coeff ( v : in out Vector; a : in coefftp ) is
  begin
    for i in v'range loop
      Mult_Coeff(v(i),a);
    end loop;
  end Mult_Coeff;

  procedure Mult_Vector ( v1 : in out Vector; v2 : in Vector ) is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then raise Range_Error;
     else for i in v1'range loop
            Mult_Coeff(v1(i),v2(i));
          end loop;
    end if;
  end Mult_Vector;

  procedure Clear ( v : in out Link_to_Vector ) is

    procedure free is new unchecked_deallocation(Vector,Link_to_Vector);

  begin
    if v /= null
     then Clear(v.all);
          free(v);
    end if;
  end Clear;

  function Equal ( v1,v2 : Link_to_Vector ) return boolean is
  begin
    if (v1 = null) and (v2 = null)
     then return true;
     elsif (v1 = null) or (v2 = null)
         then return false;
         else return Equal(v1.all,v2.all);
    end if;
  end Equal;

  procedure Copy ( v1: in Link_to_Vector; v2 : in out Link_to_Vector ) is
  begin
    Clear(v2);
    if v1 /= null
     then v2 := new Vector(v1'range);
          for i in v1'range loop
            v2(i) := v1(i);
          end loop;
    end if;
  end Copy;

  function "+" ( v1,v2 : Link_to_Vector ) return Link_to_Vector is
  begin
    if v1 = null
     then return v2;
     elsif v2 = null
         then return v1;
         else return new Vector'(v1.all + v2.all);
    end if;
  end "+";

  function "-" ( v : Link_to_Vector ) return Link_to_Vector is
  begin
    if v = null
     then return v;
     else return new Vector'(-v.all);
    end if;
  end "-";

  function "-" ( v1,v2 : Link_to_Vector ) return Link_to_Vector is
  begin
    if v2 = null
     then return v1;
     elsif v1 = null
         then return -v2;
         else return new Vector'(v1.all - v2.all);
    end if;
  end "-";

  function "*" ( v : Link_to_Vector; a : coefftp ) return Link_to_Vector is
  begin
    if v = null 
     then return null;
     else return new Vector'(v.all*a);
    end if;
  end "*";

  function "*" ( a : coefftp; v : Link_to_Vector ) return Link_to_Vector is
  begin
    return v*a;
  end "*";

  function "*" ( v1,v2 : Link_to_Vector ) return coefftp is
  begin
    if (v1 = null) or (v2 = null)
     then return zero;
     else return v1.all*v2.all;
    end if;
  end "*";

  function "*" ( v1,v2 : Link_to_Vector ) return Link_to_Vector is
  begin
    if (v1 = null) or (v2 = null)
     then return null;
     else return new Vector'(v1.all*v2.all);
    end if;
  end "*";

  function Sum ( v : Link_to_Vector ) return coefftp is
  begin
   if v = null
    then return zero;
    else return Sum(v.all);
   end if;
  end Sum;

  procedure Plus_Vector ( v1 : in out Link_to_Vector;
                          v2 : in Link_to_Vector ) is
  begin
    if v2 = null
     then null;
     elsif v1 = null
         then Copy(v2,v1);
         else Plus_Vector(v1.all,v2.all);
    end if;
  end Plus_Vector;

  procedure Min_Vector ( v : in out Link_to_Vector ) is
  begin
    if v = null
     then null;
     else Min_Vector(v.all);
    end if;
  end Min_Vector;

  procedure Min_Vector ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector ) is
  begin
    if v2 = null
     then null;
     elsif v1 = null
         then v1 := new Vector'(v2.all);
              Min_Vector(v1.all);
         else Min_Vector(v1.all,v2.all);
    end if;
  end Min_Vector;

  procedure Mult_Coeff ( v : in out Link_to_Vector; a : in coefftp ) is
  begin
    if v /= null
     then Mult_Coeff(v.all,a);
    end if;
  end Mult_Coeff;

  procedure Mult_Vector ( v1 : in out Link_to_Vector;
                          v2 : in Link_to_Vector ) is
  begin
    if v2 = null
     then null;
     elsif v1 = null
         then Clear(v1);
         else Mult_Vector(v1.all,v2.all);
    end if;
  end Mult_Vector;

end Vectors;
