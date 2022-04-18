package body Lists is

-- INTERNAL DATA :

  type Node is record
    The_Item : Item;
    Next     : List;
  end record;

  Free_List : List := null;

-- AUXILIARIES :

  procedure Set_Next ( The_Node : in out Node; To_Next : in List ) is
  begin
    The_Node.Next := To_Next;
  end Set_Next;

  function Next_Of ( The_Node : in Node ) return List is
  begin
    return The_Node.Next;
  end Next_Of;

  procedure Free ( l : in out List ) is

    tmp : List;

  begin
    while l /= null loop
      tmp := l;
      l := Next_Of(l.all);
      Set_Next(tmp.all,Free_List);
      Free_List := tmp;
    end loop;
  end Free;

  function New_Item return List is

    tmp : List;

  begin
    if Free_List = null
     then return new Node;
     else tmp := Free_List;
          Free_List := Next_Of(tmp.all);
          Set_Next(tmp.all,null);
          return tmp;
    end if;
  end New_Item;

-- CONSTRUCTORS :

  procedure Construct ( i : in Item; l : in out List ) is

    tmp : List;

  begin
    tmp := New_Item;
    tmp.The_Item := i;
    tmp.Next := l;
    l := tmp;
  exception
    when Storage_Error => raise Overflow;
  end Construct;

  procedure Append ( first,last : in out List; i : in Item ) is
  begin
    if Is_Null(first)
     then Construct(i,first);
          last := first;
     else declare
            tmp : List;
          begin
            Construct(i,tmp);
            Swap_Tail(last,tmp);
            last := Tail_Of(last);
          end;
    end if;
  end Append;

  procedure Concat ( first,last : in out List; l : in List ) is

    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      Append(first,last,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end Concat;

  procedure Set_Head ( l : in out List; i : in Item ) is
  begin
    l.The_Item := i;
  exception
    when Constraint_Error => raise List_Is_Null;
  end Set_Head;

  procedure Swap_Tail ( l1,l2 : in out List ) is

    tmp : List;

  begin
    tmp := l1.Next;
    l1.Next := l2;
    l2 := tmp;
  exception
    when Constraint_Error => raise List_Is_Null;
  end Swap_Tail;

  procedure Copy ( l1 : in List; l2 : in out List ) is

    From_Index : List := l1;
    To_Index   : List;

  begin
    Free(l2);
    if l1 /= null
     then l2 := New_Item;
          l2.The_Item := From_Index.The_Item;
          To_Index := l2;
          From_Index := From_Index.Next;
          while From_Index /= null loop
            To_Index.Next := New_Item;
            To_Index := To_Index.Next;
            To_Index.The_Item := From_Index.The_Item;
            From_Index := From_Index.Next;
          end loop;
    end if;
  exception
    when Storage_Error => raise Overflow;
  end Copy;

-- SELECTORS :

  function Is_Equal ( l1,l2 : List ) return boolean is

    left_index  : List := l1;
    right_index : List := l2;

  begin
    while left_index /= null loop
      if left_index.The_Item /= right_index.The_Item
       then return False;
      end if;
      left_index := left_index.Next;
      right_index := right_index.Next;
    end loop;
    return (right_index = null);
  exception
    when Constraint_Error => return false;
  end Is_Equal;

  function Length_Of ( l : List ) return natural is

    cnt : natural := 0;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      tmp := Tail_Of(tmp);
    end loop;
    return cnt;
  end Length_Of;
      
  function Is_Null ( l : list ) return boolean is
  begin
    return (l = null);
  end Is_Null;

  function Head_Of ( l : List ) return Item is
  begin
    return l.The_Item;
  exception
    when Constraint_Error => raise List_Is_Null;
  end Head_Of;

  function Tail_Of ( l : List ) return List is
  begin
    return l.Next;
  exception
    when Constraint_Error => raise List_Is_Null;
  end Tail_Of;

-- DESTRUCTOR :

  procedure Clear ( l : in out List ) is
  begin
    Free(l);
  end Clear;

end Lists;
