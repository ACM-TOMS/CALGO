with integer_io;                       use integer_io;
with Integer_Vectors_of_Vectors;       use Integer_Vectors_of_Vectors;
with Integer_Vectors_of_Vectors_io;    use Integer_Vectors_of_Vectors_io;

package body Integer_Faces_of_Polytope_io is

-- DESCRIPTION :
--   This package contains routines for output of faces of
--   convex polytopes.

  procedure put ( f : in Face ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Face ) is
  begin
    put(file," spanned by "); put(file,f.all'length,1);
    put_line(file," points :"); put(file,f.all);
  end put;

  procedure put ( f : in Faces ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Faces ) is

    cnt : natural := 0;
    tmp : Faces := f;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      put(file,"Face "); put(file,cnt,1); put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( f : in Array_of_Faces ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Array_of_Faces ) is
  begin
    for i in f'range loop
      if not Is_Null(f(i))
       then put(file,"faces at component "); put(file,i,1);
            put_line(file," :"); put(file,f(i));
      end if;
    end loop;
  end put;

end Integer_Faces_of_Polytope_io;
