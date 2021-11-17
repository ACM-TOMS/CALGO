with integer_io;    use integer_io;

package body Integer_Matrices_io is

  procedure get ( m : out Matrix ) is
  begin
    get(Standard_Input,m);
  end get;

  procedure get ( m : out Matrix; rw1,rw2 : in integer ) is
  begin
    get(Standard_Input,m,rw1,rw2);
  end get;

  procedure get ( file : in file_type; m : out Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        get(file,m(i,j));
      end loop;
    end loop;
  end get;

  procedure get ( file : in file_type;
                  m : out Matrix; rw1,rw2 : in integer ) is
  begin
    for i in rw1..rw2 loop
      for j in m'range(2) loop
        get(file,m(i,j));
      end loop;
    end loop;
  end get;

  procedure put ( m : in Matrix ) is
  begin
    put(Standard_Output,m);
  end put;

  procedure put ( m : in Matrix; width : in natural ) is
  begin
    put(Standard_Output,m,width);
  end put;

  procedure put ( m : in Matrix; rw1,rw2 : in integer ) is
  begin
    put(Standard_Output,m,rw1,rw2);
  end put;

  procedure put ( file : in file_type; m : in Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j),1);
      end loop;
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type; m : in Matrix; width : in natural ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j),width);
      end loop;
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type; m : in Matrix; rw1,rw2 : in integer ) is
  begin
    for i in rw1..rw2 loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j),1);
      end loop;
      new_line(file);
    end loop;
  end put;

end Integer_Matrices_io;
