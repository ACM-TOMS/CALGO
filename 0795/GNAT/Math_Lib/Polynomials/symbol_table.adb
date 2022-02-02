with unchecked_deallocation;

package body Symbol_Table is

-- INTERNAL DATA :

  type Symbol_Array is array(positive range <>) of Symbol;
  type Symbol_Table ( max : natural ) is record
    number : natural;     -- number of symbols that are not blank
    syms : Symbol_Array(1..max);
  end record;
  type Link_to_Symbol_Table is access Symbol_Table;

  st : Link_to_Symbol_Table;

-- CREATORS :

  procedure Init ( max : in natural ) is
  begin
    st := new Symbol_Table(max);
    st.all.number := 0;
  end Init;

  procedure Enlarge ( max : in natural ) is
  begin
    if Empty
     then Init(max);
     else declare
            oldst : Symbol_Array(1..st.number);
            maxst : constant natural := max + st.max;
          begin
            for i in oldst'range loop
              oldst(i) := st.syms(i);
            end loop;
            Clear;
            Init(maxst);
            for i in oldst'range loop
              Add(oldst(i));
            end loop;
          end;
    end if;
  end Enlarge;

-- CONSTRUCTOR :

  procedure Add ( sb : in Symbol; pos : out natural ) is

    tab : symbol_table renames st.all;

  begin
    tab.number := tab.number + 1;
    for i in sb'range loop
      tab.syms(tab.number)(i) := sb(i);
    end loop;
    pos := tab.number;
  exception
    when others => raise OVERFLOW_IN_THE_SYMBOL_TABLE;
  end Add;

  procedure Add ( sb : in Symbol ) is

    pos : natural;

  begin
    Add(sb,pos);
  end Add;

-- SELECTORS :

  function "<" ( s1,s2 : Symbol ) return boolean is
  begin
    for i in s1'range loop
      if s1(i) < s2(i)
       then return true;
       elsif s1(i) > s2(i)
           then return false;
      end if;
    end loop;
    return false;
  end "<";

  function ">" ( s1,s2 : Symbol ) return boolean is
  begin
    for i in s1'range loop
      if s1(i) > s2(i)
       then return true;
       elsif s1(i) < s2(i)
           then return false;
      end if;
    end loop;
    return false;
  end ">";

  function Number return natural is
  begin
    if st = null
     then return 0;
     else return st.all.number;
    end if;
  end Number;
  
  function Empty return boolean is
  begin
    return (st = null);
  end Empty;

  function Get ( sb : symbol ) return natural is

    tab : symbol_table renames st.all;

  begin
    for i in 1..tab.number loop
      if tab.syms(i) = sb
       then return i;
      end if;
    end loop;
    return 0;
  end Get;

  function Get ( i : natural ) return symbol is

    tab : symbol_table renames st.all;

  begin
    if i > tab.number
     then raise INDEX_OUT_OF_RANGE;
     else return tab.syms(i);
    end if;
  end Get;

-- DESTRUCTOR :

  procedure Clear is

    procedure free is 
      new unchecked_deallocation (symbol_table,link_to_symbol_table);

  begin
    free(st);
  end Clear;

end Symbol_Table;
