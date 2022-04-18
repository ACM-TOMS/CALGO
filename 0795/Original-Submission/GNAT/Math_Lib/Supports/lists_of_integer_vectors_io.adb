with Integer_Vectors;          use Integer_Vectors;
with Integer_Vectors_io;       use Integer_Vectors_io;

package body Lists_of_Integer_Vectors_io is

  procedure get ( n,m : in natural; l : out List ) is
  begin
    get(Standard_Input,n,m,l);
  end get;

  procedure get ( file : in file_type; n,m : in natural; l : out List ) is

    tmp,tmp_last : List;

  begin
    for i in 1..m loop
      declare
	v : Link_to_Vector;
      begin
        get(file,n,v);
        Append(tmp,tmp_last,v);
      end;
    end loop;
    l := tmp;
  end get;

  procedure put ( l : in List ) is
  begin
    put(Standard_Output,l);
  end put;

  procedure put ( file : in file_type; l : in List ) is

    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end put;

end Lists_of_Integer_Vectors_io;
