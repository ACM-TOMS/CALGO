with integer_io;                     use integer_io;
with Communications_with_User;       use Communications_with_User;
with Symbol_Table,Symbol_Table_io;   use Symbol_Table;

package body Symbolic_Symmetry_Group_io is

  procedure iget ( p : in out Permutation ) is

  -- DESCRIPTION ;
  --   Auxiliary procedure to realize an interactive get.

    ok : boolean := false;
    pp : Permutation(p'range) := (p'range => 0);
    ans : character;

  begin
    while not ok loop
      get(Standard_Input,pp);
      ok := Is_Permutation(pp);
      if not ok
       then put("is no permutation.  Do you want to retry ? (y/n) ");
            Ask_Yes_or_No(ans);
            ok := (ans /= 'y');
      end if;
    end loop;
    p := pp;
  end iget;

  procedure get ( p : in out Permutation ) is
  begin
    iget(p);
  end get;

  procedure get ( file : in file_type; p : in out Permutation ) is
 
    sb : Symbol;
    ch : character;

  begin
    for i in p'range loop
      sb := (sb'range => ' ');
      get(ch);
      while ch = ' ' loop get(ch); end loop;
      if ch = '-'
       then get(ch);
            Symbol_Table_io.get(file,ch,sb);
            p(i) := -Symbol_Table.get(sb);
       else Symbol_Table_io.get(file,ch,sb);
            p(i) := Symbol_Table.get(sb);
      end if;
    end loop;
  end get;

  procedure get ( l : in out List_of_Permutations; n,nb : in natural ) is

    p : Permutation(1..n) := (1..n => 0);
    l2 : List_of_Permutations;

  begin
    for i in 1..nb loop
      iget(p);
      Append(l,l2,p);
    end loop;
  end get;

  procedure get ( file : in file_type;
		  l : in out List_of_Permutations; n,nb : in natural ) is

    p : Permutation(1..n) := (1..n => 0);
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

    sb : Symbol;

  begin
    for i in p'range loop
      put(file,' ');
      if p(i) < 0
       then put(file,'-');
            sb := Symbol_Table.get(-p(i));
       else sb := Symbol_Table.get(p(i));
      end if;
      Symbol_Table_io.put(file,sb);
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

end Symbolic_Symmetry_Group_io;
