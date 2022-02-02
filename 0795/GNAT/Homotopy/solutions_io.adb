with integer_io;                     use integer_io;
with Floating_Point_Numbers;
with Communications_with_User;       use Communications_with_User;
with Complex_Numbers;                use Complex_Numbers;
with Complex_Numbers_io;             use Complex_Numbers_io;
with Complex_Vectors_io;             use Complex_Vectors_io;
with Symbol_Table;                   use Symbol_Table;

package body Solutions_io is

  use Floating_Point_Numbers.double_float_io;

-- INPUT OF SYMBOL :

  procedure skip_symbol ( file : in file_type ) is

  -- DESCRIPTION :
  --   Skips all symbols until a `:' is encountered.

    c : character;

  begin
    loop
      get(file,c);
      exit when (c = ':');
    end loop;
  end skip_symbol;

  function get_symbol ( file : in file_type ) return natural is

  -- DESCRIPTION :
  --   Reads a symbol from standard input and returns its number.

    sb : Symbol := (1..3 => ' ');
    c : character;

  begin
    loop       -- skip the spaces
      get(file,c);
      exit when c /= ' ';
    end loop;
    sb(1) := c;
    for i in 2..3 loop
      get(file,c);
      exit when c = ' ';
      sb(i) := c;
    end loop;
    return Symbol_Table.get(sb);
  end get_symbol;    

-- OUTPUT OF A SYMBOL :

  procedure put_symbol ( file : in file_type; i : in natural ) is

  -- DESCRIPTION :
  --   Given the number of the symbol,
  --   the corresponding symbol will be written.

    sb : Symbol := Get(i);

  begin
    for k in sb'range loop
      exit when sb(k) = ' ';
      put(file,sb(k));
    end loop;
  end put_symbol;

-- INPUT OF A SOLUTION VECTOR :

  procedure get_vector ( s : in out Solution ) is
  begin
    get_vector(Standard_Input,s);
  end get_vector;

  procedure get_vector ( file : in file_type; s : in out Solution ) is

    ind : natural;

  begin
    if Symbol_Table.Number < s.n
     then for i in s.v'range loop
            skip_symbol(file); get(file,s.v(i));
          end loop;
     else for i in s.v'range loop
            ind := get_symbol(file); skip_symbol(file);
            get(file,s.v(ind));
          end loop;
    end if; 
  end get_vector;

-- INPUT OF A SOLUTION :

  procedure get ( s : in out Solution ) is
  begin
    get(Standard_Input,s);
  end get;

  procedure get ( file : in file_type; s : in out Solution ) is

    c : character;

  begin
    get(file,c); get(file,c); get(file,c); get(file,c); get(file,s.t);
    get(file,c); get(file,c); get(file,c); get(file,c); get(file,s.m);
    if not End_of_Line(file)
     then get(file,c); Skip_line(file);  -- skip information on this line
    end if;
    get(file,c); skip_line(file);
    get_vector(file,s);
  end get;

-- OUTPUT OF A SOLUTION VECTOR :

  procedure put_vector ( s : in Solution ) is
  begin
    put_vector(Standard_Output,s);
  end put_vector;

  procedure put_vector ( file : in file_type; s : in Solution ) is
  begin
    if Symbol_Table.Number < s.n
     then for i in s.v'range loop
            put(file," x"); put(file,i,1); put(file," : ");
            put(file,s.v(i)); new_line(file);
          end loop;
     else for i in s.v'range loop
            put(file,' '); put_symbol(file,i); put(file," : ");
            put(file,s.v(i)); new_line(file);
          end loop;
    end if;
  end put_vector;

-- OUTPUT OF A SOLUTION :

  procedure put ( s : in Solution ) is
  begin
    put(Standard_Output,s);
  end put;

  procedure put ( file : in file_type; s : in Solution ) is
  begin
    put(file,"t : "); put(file,s.t); new_line(file);
    put(file,"m : "); put(file,s.m,1); new_line(file);
    put_line(file,"the solution for t :");
    put_vector(file,s);
    put(file,"==");
    put(file," err : "); put(file,s.err,2,3,3); put(file," =");
    put(file," rco : "); put(file,s.rco,2,3,3); put(file," =");
    put(file," res : "); put(file,s.res,2,3,3); put(file," =");
  end put;

-- INPUT OF A LIST OF SOLUTIONS :

  procedure get ( len,n : in natural;
                  sols,sols_last : in out Solution_List ) is
  begin
    get(Standard_Input,sols,sols_last);
  end get;

  procedure get ( len,n : in natural; sols : in out Solution_List ) is
  begin
    get(Standard_Input,len,n,sols);
  end get;

  procedure get ( sols,sols_last : in out Solution_List ) is
  begin
    get(Standard_Input,sols,sols_last);
  end get;

  procedure get ( sols : in out Solution_List ) is
  begin
    get(Standard_Input,sols);
  end get;

  procedure get ( file : in file_type; len,n : in natural;
                  sols,sols_last : in out Solution_List ) is

    s : Solution(n);
    c : character;

  begin
    for i in 1..len loop
      get(file,c); skip_line(file);    -- skip opening bar
      get(file,c); skip_line(file);    -- skip line with solution number
      get(file,s);
      Append(sols,sols_last,s);
    end loop;
    get(file,c); skip_line(file);     -- skip closing bar
  end get;

  procedure get ( file : in file_type; len,n : in natural;
                  sols : in out Solution_List ) is

    sols_last : Solution_List;

  begin
    get(file,len,n,sols,sols_last);
  end get;

  procedure get ( file : in file_type;
                  sols,sols_last : in out Solution_List ) is

    len,n : natural;

  begin
    get(file,len); get(file,n);
    get(file,len,n,sols,sols_last);
  end get;

  procedure get ( file : in file_type; sols : in out Solution_List ) is

    len,n : natural;

  begin
    get(file,len); get(file,n);
    get(file,len,n,sols);
  end get;

-- OUTPUT OF A LIST OF SOLUTIONS :

  procedure put_bar ( file : in file_type ) is
  begin
    put_line(file,
             "===========================================================");
  end put_bar;

  procedure put ( sols : in Solution_List ) is
  begin
    put(Standard_Output,sols);
  end put;

  procedure put ( len,n : in natural; sols : in Solution_List ) is
  begin
    put(Standard_Output,len,n,sols);
  end put;

  procedure put ( file : in file_type; sols : in Solution_List ) is
  begin
    if not Is_Null(sols)
     then declare
            count : natural := 1;
            temp : Solution_List := sols;
          begin
            put_bar(file);
            while not Is_Null(temp) loop
              put(file,"solution "); put(file,count,1);
              put(file," :"); new_line(file);
              put(file,Head_Of(temp).all);
              put_line(file,"="); -- instead of : put_bar(file);
              temp := Tail_Of(temp);
              count := count + 1;
            end loop;
          end;
    end if;
  end put;

  procedure put ( file : in file_type; len,n : in natural;
                  sols : in Solution_List ) is
  begin
    put(file,len,1); put(file," "); put(file,n,1); new_line(file);
    put(file,sols);
  end put;

  procedure Display_Format is

    s : array(1..24) of string(1..65);

  begin
    s( 1):="  A solution list of a complex polynomial system  is  denoted  by";
    s( 2):="the  number of solutions and the dimension, followed by a list of";
    s( 3):="solutions.   The  solutions  are  separated  by  a  banner  line,";
    s( 4):="followed by their position in the list.                          ";
    s( 5):="  A solution consists of the current value  of  the  continuation";
    s( 6):="parameter  t,  its  multiplicity  (or  winding number) m, and the";
    s( 7):="solution vector.                                                 ";
    s( 8):="  A solution vector contains as many lines as the dimension.  The";
    s( 9):="i-th  line  starts  with  the  symbol  that  represents  the i-th";
    s(10):="unknown, followed by the colon `:' and two floating-point numbers";
    s(11):="representing  respectively  the  real  and  imaginary part of the";
    s(12):="solution component.                                              ";
    s(13):="  As example we list the solution  list of  the  regular solution";
    s(14):="(1,i) of a 2-dimensional system in the unknowns x and y at t=1.  ";
    s(15):="                                                                 ";
    s(16):="1 2                                                              ";
    s(17):="=================================================================";
    s(18):="solution 1 :                                                     ";
    s(19):="t :  1.00000000000000E+00  0.00000000000000E+00                  ";
    s(20):="m : 1                                                            ";
    s(21):="the solution for t :                                             ";
    s(22):=" x : 1.00000000000000E+00  0.00000000000000E+00                  ";
    s(23):=" y : 0.00000000000000E+00  1.00000000000000E+00                  ";
    s(24):="=================================================================";
    for i in s'range loop
      put_line(s(i));
    end loop;
  end Display_Format;

  procedure Read ( sols : in out Solution_List ) is

    file : file_type;

  begin
    put_line("Reading the name of the file for the solutions.");
    Read_Name_and_Open_File(file);
    get(file,sols);
    Close(file);
  exception
    when others => Close(file); Clear(sols);
                   put_line("INCORRECT FORMAT OF SOLUTION LIST");
                   Display_Format; new_line;
                   Read(sols);
  end Read;

end Solutions_io;
