with integer_io;   use integer_io;

package body Symmetry_Group_io is

  procedure get ( p : out Permutation ) is
  begin
    get(Standard_Input,p);
  end get;

  procedure get ( file : in file_type; p : out Permutation ) is
  begin
    for i in p'range loop
      get(file,p(i));
    end loop;
  end get;

  procedure get ( l : in out List_of_Permutations; n,nb : in natural ) is
  begin
    get(Standard_Input,l,n,nb);
  end get;

  procedure get ( file : in file_type;
		  l : in out List_of_Permutations; n,nb : in natural ) is

    p : Permutation(1..n);
    l2 : List_of_Permutations;

  begin
    for i in 1..nb loop
      get(file,p);
      if Is_Permutation(p)
       then Append(l,l2,p);
      end if;
    end loop;
  end get;
      
  procedure put ( p : in Permutation ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( file : in file_type; p : in Permutation ) is
  begin
    for i in p'range loop
      put(file,' '); put(file,p(i),1);
    end loop;
  end put;

  procedure put ( l : in List_of_Permutations ) is
  begin
    put(Standard_Output,l);
  end put;

  procedure put ( file : in file_type; l : in List_of_Permutations ) is

    temp : List_of_Permutations := l;

  begin
    while not Is_Null(temp) loop
      put(file,Permutation(Head_Of(temp).all));
      new_line(file);
      temp := Tail_Of(temp);
    end loop;
  end put;

end Symmetry_Group_io;
