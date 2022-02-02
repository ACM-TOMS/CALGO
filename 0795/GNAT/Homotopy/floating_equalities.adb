package body Floating_Equalities is

  function Is_Equal ( x,y,tol : double_float ) return boolean is
  begin
    if abs(x-y) <= tol
     then return true;
     else return false;
    end if;
  end Is_Equal;

  function Is_Equal ( x,y : double_complex; tol : double_float )
                    return boolean is
  begin
    if Is_Equal(REAL_PART(x),REAL_PART(y),tol)
      and Is_Equal(IMAG_PART(x),IMAG_PART(y),tol)
     then return true;
     else return false;
    end if;
  end Is_Equal;

  function Is_Equal ( x,y : Float_Vectors.Vector; tol : double_float )
                    return boolean is
  begin
    for i in x'range loop
      if not Is_Equal(x(i),y(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Is_Equal ( x,y : Complex_Vectors.Vector; tol : double_float )
                    return boolean is
  begin
    for i in x'range loop
      if not Is_Equal(x(i),y(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

end Floating_Equalities;
