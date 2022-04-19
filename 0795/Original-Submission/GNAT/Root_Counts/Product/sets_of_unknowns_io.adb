with integer_io;                           use integer_io;
with Symbol_Table,Symbol_Table_io;         use Symbol_Table;

package body Sets_of_Unknowns_io is

  procedure get ( s : in out Set ) is
  begin
    get(Standard_Input,s);
  end get;

  procedure get ( file : in file_type; s : in out Set ) is

    ch : character;
    continue : boolean;

    procedure Read_Symbol ( cont : out boolean ) is
    begin
      loop
	get(file,ch);
	exit when ch /= ' ';
      end loop;
      if ch = '}'
       then cont := false;
       else declare
	      sb : symbol;
	      nb : natural;
            begin
              sb := (sb'range => ' ');
              Symbol_Table_io.get(file,ch,sb,'}');
	      nb := Symbol_Table.Get(sb);
	      if nb = 0
	       then Symbol_Table.Add(sb,nb);
              end if;
	      Add(s,nb);
            end;
	    cont := (ch /= '}');
      end if;
    end Read_Symbol;

  begin
    loop
      get(file,ch);
      exit when ch = '{' or ch = '}';
    end loop;
    if ch = '{'
     then loop
            Read_Symbol(continue);
            exit when not continue;
          end loop;
    end if;
  end get;

  procedure put ( s : in Set ) is
  begin
    put(Standard_Output,s);
  end put;

  procedure put ( file : in file_type; s : in Set ) is
 
    standard : boolean := (Symbol_Table.Number = 0);
 
    procedure Write_Element ( i: in natural ) is
    begin
      if standard
       then put(file,'x'); put(file,i,2);
            if i >= 10
             then put(file,' ');
            end if;
       else Symbol_Table_io.put(file,Symbol_Table.Get(i));
            put(file,' ');
      end if;
    end Write_Element;

  begin
    text_io.put(file,"{");
    for i in 1..Dimension(s) loop
      if Is_In(s,i)
       then Write_Element(i);
      end if;
    end loop;
    text_io.put(file,"}");
  end put;

end Sets_of_Unknowns_io;
