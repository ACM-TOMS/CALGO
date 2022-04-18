with Integer_Vectors_io;    use Integer_Vectors_io;

package body Integer_Vectors_of_Vectors_io is

  procedure get ( m : in natural; v : in out Vector ) is
  begin
    get(Standard_Input,m,v);
  end get;

  procedure get ( file : in file_type; m : in natural; v : in out Vector ) is
  begin
    for i in v'range loop
      get(file,m,v(i));
    end loop;
  end get;

  procedure get ( n,m : in natural; v : out Link_to_Vector ) is
  begin
    get(Standard_Input,n,m,v);
  end get;

  procedure get ( file : in file_type; n,m : in natural;
                  v : out Link_to_Vector ) is
    vv : Vector(1..n);
  begin
    for i in vv'range loop
      get(file,m,vv(i));
    end loop;
    v := new Vector'(vv);
  end get;

  procedure put ( v : in Vector ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Vector ) is
  begin
    for i in v'range loop
      put(file,v(i)); new_line(file);
    end loop;
  end put;

  procedure put ( v : in Link_to_Vector ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_Vector ) is
  begin
    if v /= null
     then put(file,v.all);
    end if;
  end put;

end Integer_Vectors_of_Vectors_io;
