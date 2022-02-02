with integer_io,Integer_Vectors_io;   use integer_io,Integer_Vectors_io;
with Lists_of_Float_Vectors_io;       use Lists_of_Float_Vectors_io;

package body Arrays_of_Float_Vector_Lists_io is

  procedure get ( al : in out Link_to_Array_of_Lists ) is
  begin
    get(Standard_Input,al);
  end get;

  procedure get ( n : in natural; m : in Vector; al : out Array_of_Lists ) is
  begin
    get(Standard_Input,n,m,al);
  end get;

  procedure get ( file : in file_type; al : in out Link_to_Array_of_Lists ) is

    n : natural;

  begin
    get(file,n);
    al := new Array_of_Lists(1..n);
    declare
      m : Vector(1..n) := (1..n => 0);
    begin
      get(file,m);
      get(file,n,m,al.all);
    end;
  end get;

  procedure get ( file : in file_type; n : in natural; m : in Vector;
		  al : out Array_of_Lists ) is
  begin
    for i in al'range loop
      get(file,n,m(i),al(i));
    end loop;
  end get;

  procedure put ( al : in Array_of_Lists ) is
  begin
    put(Standard_Output,al);
  end put;

  procedure put ( file : in file_type; al : in Array_of_Lists ) is
  begin
    for i in al'range loop
      put(file,al(i)); new_line(file);
    end loop;
  end put;

  procedure put ( al : in Array_of_Lists; fore,aft,exp : in natural ) is
  begin
    put(Standard_Output,al,fore,aft,exp);
  end put;

  procedure put ( file : in file_type;
                  al : in Array_of_Lists; fore,aft,exp : in natural ) is
  begin
    for i in al'range loop
      put(file,al(i),fore,aft,exp); new_line(file);
    end loop;
  end put;

end Arrays_of_Float_Vector_Lists_io;
