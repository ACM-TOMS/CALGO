with integer_io,Numbers_io;             use integer_io,Numbers_io;
with Communications_with_User;          use Communications_with_User;
with Symbol_Table,Symbol_Table_io;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package body Complex_Polynomial_Systems_io is

-- SCANNING THE LINE FOR A NATURAL NUMBER :

  function Scan_Line ( file : in file_type ) return natural is

    m : natural := 0;
    ch : character;

  begin
    while not END_OF_LINE(file) loop
      get(file,ch);
      case ch is
        when '0' => m := 10*m;
        when '1' => m := 10*m + 1;
        when '2' => m := 10*m + 2;
        when '3' => m := 10*m + 3;
        when '4' => m := 10*m + 4;
        when '5' => m := 10*m + 5;
        when '6' => m := 10*m + 6;
        when '7' => m := 10*m + 7;
        when '8' => m := 10*m + 8;
        when '9' => m := 10*m + 9;
        when others => null;
      end case;
    end loop;
    return m;
  end Scan_Line;

-- EXCEPTION HANDLERS :

  procedure Write_Symbol_Table is

  -- DESCRIPTION :
  --   Writes the current list of symbols on one line on standard output.

  begin
    put("Current symbols : ");
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(i)); put(" ");
    end loop;
    new_line;
  end Write_Symbol_Table;

  procedure Handler ( k : in natural ) is
  begin
    put(" raised while reading polynomial "); put(k,1);
    put_line(".");
  end Handler;

-- THE INPUT OPERATIONS :

  procedure get ( n : in natural; s : out Poly_Sys ) is
  begin
    get(Standard_Input,n,s);
  end get;

  procedure get ( n,m : in natural; s : out Poly_Sys ) is
  begin
    get(Standard_Input,n,m,s);
  end get;

  procedure get ( file : in file_type; n : in natural; s : out Poly_Sys ) is

    i : natural := 1;

  begin
    while i <= n loop
      get(file,n,s(i));
      i := i + 1;
    end loop;
  exception
    when ILLEGAL_CHARACTER    => put("ILLEGAL_CHARACTER");    Handler(i); raise;
    when ILLEGAL_SYMBOL       => put("ILLEGAL_SYMBOL");       Handler(i); raise;
    when ILLEGAL_OPERATION    => put("ILLEGAL_OPERATION");    Handler(i); raise;
    when INFINITE_NUMBER      => put("INFINITE_NUMBER");      Handler(i); raise;
    when OVERFLOW_OF_UNKNOWNS => put("OVERFLOW_OF_UNKNOWNS"); Handler(i);
                                 Write_Symbol_Table; raise;
    when BAD_BRACKET          => put("BAD_BRACKET");          Handler(i); raise;
  --when others               => put("Exception");            Handler(i); raise;
  end get;

  procedure get ( file : in file_type; n,m : in natural; s : out Poly_Sys ) is

    i : natural := 1;

  begin
    while i <= n loop
      get(file,m,s(i));
      i := i + 1;
    end loop;
  exception
    when ILLEGAL_CHARACTER    => put("ILLEGAL_CHARACTER");    Handler(i); raise;
    when ILLEGAL_SYMBOL       => put("ILLEGAL_SYMBOL");       Handler(i); raise;
    when ILLEGAL_OPERATION    => put("ILLEGAL_OPERATION");    Handler(i); raise;
    when INFINITE_NUMBER      => put("INFINITE_NUMBER");      Handler(i); raise;
    when OVERFLOW_OF_UNKNOWNS => put("OVERFLOW_OF_UNKNOWNS"); Handler(i);
                                 Write_Symbol_Table; raise;
    when BAD_BRACKET          => put("BAD_BRACKET");          Handler(i); raise;
  --when others               => put("Exception");            Handler(i); raise;
  end get;

  procedure get ( s : out Poly_Sys ) is
  begin
    get(Standard_Input,s);
  end get;

  procedure get ( file : in file_type; s : out Poly_Sys ) is

    n,m : natural;

  begin
    get(file,n);
    m := Scan_Line(file);
    if m /= 0
     then get(file,n,m,s);
     else get(file,n,s);
    end if;
  end get;

-- MORE USER FRIENDLY INPUT OPERATIONS :

  procedure get ( lp : in out Link_to_Poly_Sys ) is

    inpt : file_type;
    n,m : natural;
    onfile : character;

  begin
  --------------------------------------------
  --  GETTING THE DIMENSION OF THE PROBLEM  --
  --------------------------------------------
    loop
      put("Is the system on a file ? (y/n/i=info) ");
      Ask_Alternative(onfile,"yni");
      if onfile = 'i'
       then new_line; Display_Format; new_line;
      end if;
      exit when onfile /= 'i';
    end loop;
    new_line;
    if onfile = 'y'
     then declare
            procedure Read is
            begin	 
              put_line("Reading the name of the input file.");
              Read_Name_and_Open_File(inpt);
              get(inpt,n);
            end Read;
          begin
	    Read;
          exception
            when others => put_line("The data on the file is not correct.");
                           put_line("A natural number is expected first.");
	                   put_line("Supply another file name."); Close(inpt);
                           Read;
          end;
     else put("Give the number of polynomials : "); Read_Natural(n);
    end if;
  ---------------------------------------
  --  GETTING THE POLYNOMIAL SYSTEM :  --
  ---------------------------------------
    lp := new Poly_Sys(1..n);
    declare
      procedure Read is
      begin 
        if onfile = 'y'
         then m := Scan_Line(inpt);
              if m /= 0
               then get(inpt,n,m,lp.all);
               else get(inpt,n,lp.all);
              end if;
              Close(inpt);
         else put("Give the number of unknowns : "); Read_Natural(m);
              put("Give "); put(n,2);
              if n = 1
               then put_line(" polynomial : ");
               else put_line(" polynomials : ");
              end if;
              if m = n
               then get(n,lp.all);
               else get(n,m,lp.all);
              end if;
              skip_line;  -- skip end_of_line symbol
        end if;
      exception
        when others => if onfile = 'y' then Close(inpt); end if;
                       put_line("Polynomial system read : "); put(lp.all,'*');
                       raise;
      end Read;
    begin
      Read;
    exception
      when others =>
             if onfile = 'y'
              then put_line("The polynomials on the file are incorrect."
                          & " Try again...");
              else put_line("The polynomials are incorrect. Try again...");
             end if;
             Clear(lp); Symbol_Table.Clear;
             get(lp);
    end;

  end get;

  procedure get ( file : in file_type; lp : in out Link_to_Poly_Sys ) is

    n,m : natural;

  begin
    get(file,n);
    lp := new Poly_Sys(1..n);
    m := Scan_Line(file);
    if m /= 0
     then get(file,n,m,lp.all);
     else get(file,n,lp.all);
    end if;
    skip_line(file);    -- skip end_of_line symbol
  end get;

-- THE OUTPUT OPERATIONS :

  procedure put ( n : out natural; s : in Poly_Sys; pow : in power := '*' ) is
  begin
    put(Standard_Output,n,s,pow);
  end put;

  procedure put ( n,m : out natural; s : in Poly_Sys; pow : in power := '*' ) is
  begin
    put(Standard_Output,n,m,s,pow);
  end put;
 
  procedure put ( file : in file_type; n : out natural; s : in Poly_Sys;
                  pow : in power := '*' ) is

    nn : natural := s'length;

  begin
    n := nn;
    for i in s'range loop
      put(file,nn,s(i),pow);
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type; n,m : out natural; s : in Poly_Sys;
                  pow : in power := '*' ) is

    nn : natural := s'length;
    mm : natural := Number_of_Unknowns(s(s'first));

  begin
    n := nn;
    m := mm;
    for i in s'range loop
      put(file,mm,s(i),pow);
      new_line(file);
    end loop;
  end put;
  
  procedure put ( s : in Poly_Sys; pow : in power ) is
  begin
    put(Standard_Output,s,pow);
  end put;

  procedure put ( file : in file_type; s : in Poly_Sys;
                  pow : in power ) is

    n : natural := s'length;
    m : natural := Number_of_Unknowns(s(s'first));

  begin
    put(file,n,2);
    if m /= n
     then put(file,' '); put(file,m,2);
    end if;
    new_line(file);
    put(file,n,s,pow);
  end put;

  procedure put ( s : in Poly_Sys ) is
  begin
    put(Standard_Output,s,'*');
  end put;

  procedure put ( file : in file_type; s : in Poly_Sys ) is
  begin
    put(file,s,'*');
  end put;

  procedure put_line ( s : in Poly_Sys ) is
  begin
    put_line(Standard_Output,s);
  end put_line;

  procedure put_line ( file : in file_type; s : in Poly_Sys ) is
  begin
    put_line(file,s,'*');
  end put_line;

  procedure put_line ( s : in Poly_Sys; pow : in Power ) is
  begin
    put_line(Standard_Output,s,pow);
  end put_line;

  procedure put_line ( file : in file_type; s : in Poly_Sys; pow : in Power ) is
  begin
    put(file,s'length,2); new_line(file);
    for i in s'range loop
      put_line(file,s(i),pow);
    end loop;
  end put_line;

  procedure Display_Format is

    s : array(1..3) of string(1..65);

  begin
    s(1):="  A  complex  polynomial  system  is  denoted  by  the  dimension";
    s(2):="followed  by  as  many  complex  multivariate  polynomials as the";
    s(3):="dimension.  The dimension is a positive natural number.          ";
    for i in s'range loop
      put_line(s(i));
    end loop;
    Complex_Multivariate_Polynomials_io.Display_Format;
  end Display_Format;

end Complex_Polynomial_Systems_io;
