package body Symbol_Table_io is

  procedure get ( sb : in out symbol ) is
  begin
    Symbol_Table_io.get(Standard_Input,sb);
  end get;

  procedure get ( sb : in out symbol; delimiter : in character ) is
  begin
    Symbol_Table_io.get(Standard_Input,sb,delimiter);
  end get;

  procedure get ( ch : in out character; sb : in out symbol ) is
  begin
    Symbol_Table_io.get(Standard_Input,ch,sb);
  end get;

  procedure get ( ch : in out character; sb : in out symbol;
                  delimiter : in character ) is
  begin
    Symbol_Table_io.get(Standard_Input,ch,sb,delimiter);
  end get;

  procedure get ( file : in file_type; sb : in out symbol ) is

    ch : character;

  begin
    text_io.get(file,ch);
    Symbol_Table_io.get(file,ch,sb);
  end get;

  procedure get ( file : in file_type; sb : in out symbol;
                  delimiter : in character ) is

    ch : character;

  begin
    text_io.get(file,ch);
    Symbol_Table_io.get(file,ch,sb,delimiter);
  end get;

  procedure get ( file : in file_type; ch : in out character; 
                  sb : in out symbol ) is

    cnt : natural;

  begin
    if ch = ' '
     then while ch = ' ' loop 
            text_io.get(file,ch); 
          end loop;
          if ch /= ' '
           then sb(sb'first) := ch; cnt := 2;
          end if;
     else sb(sb'first) := ch; cnt := 2;
    end if;
    if not End_of_Line(file)
     then
       for i in sb'first+1..sb'last loop
         text_io.get(file,ch); cnt := i;
         exit when ch = ' ';
         sb(i) := ch; cnt := i+1;
         exit when End_of_Line(file);
       end loop;
    end if;
    sb(cnt..sb'last) := (cnt..sb'last => ' ');
  end get;

  procedure get ( file : in file_type; ch : in out character;
                  sb : in out symbol; delimiter : in character ) is

    cnt : natural;

  begin
    if ch = ' '
     then while ch = ' ' loop
            text_io.get(file,ch);
          end loop;
          if ch /= ' '
           then sb(sb'first) := ch; cnt := 2;
          end if;
     else sb(sb'first) := ch; cnt := 2;
    end if;
    if not End_of_Line(file)
     then
       for i in sb'first+1..sb'last loop
         text_io.get(file,ch); cnt := i;
         exit when ch = ' ' or else ch = delimiter;
         sb(i) := ch; cnt := i+1;
         exit when End_of_Line(file);
       end loop;
    end if;
    sb(cnt..sb'last) := (cnt..sb'last => ' ');
  end get;

  procedure put ( sb : in symbol ) is
  begin
    for i in sb'range loop
      exit when sb(i) = ' ';
      text_io.put(sb(i));
    end loop;
  end put;

  procedure put ( file : in file_type; sb : in symbol ) is
  begin
    for i in sb'range loop
      exit when sb(i) = ' ';
      text_io.put(file,sb(i));
    end loop;
  end put;

end Symbol_Table_io;
