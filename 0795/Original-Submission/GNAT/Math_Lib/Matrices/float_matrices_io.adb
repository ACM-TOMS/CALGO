with Floating_Point_Numbers;

package body Float_Matrices_io is

  use Floating_Point_Numbers.double_float_io;

  procedure get ( m : out Matrix ) is
  begin
    get(Standard_Input,m);
  end get;

  procedure get ( file : in file_type; m : out Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        get(file,m(i,j));
      end loop;
    end loop;
  end get;

  procedure put ( m : in Matrix ) is
  begin
    put(Standard_Output,m);
  end put;

  procedure put ( file : in file_type; m : in Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j));
      end loop;
      new_line(file);
    end loop;
  end put;

  procedure put ( m : in Matrix; fore,aft,exp : in natural ) is
  begin
    put(Standard_Output,m,fore,aft,exp);
  end put;

  procedure put ( file : in file_type; m : in Matrix;
                  fore,aft,exp : in natural ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j),fore,aft,exp);
      end loop;
      new_line(file);
    end loop;
  end put;

end Float_Matrices_io;
