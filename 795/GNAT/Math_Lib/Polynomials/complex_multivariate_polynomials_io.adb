with integer_io;                        use integer_io;
with Strings_to_Natural_Numbers;        use Strings_to_Natural_Numbers;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers,Natural_Vectors;   use Complex_Numbers;
with Symbol_Table,Symbol_Table_io;      use Symbol_Table;

package body Complex_Multivariate_Polynomials_io is

  use Floating_Point_Numbers.double_float_io;

-- INTERNAL VARIABLES :

  right : constant natural := 75;    -- these variables are needed
  column : natural := 0;             -- for the output of long polynomials

  procedure init_line is
  begin
    column := 0;
  end init_line;

  procedure line ( file : in file_type; n : natural ) is

  -- DESCRIPTION :
  --   this routine decides when a new line on the output has to be taken;
  --   n is the number of characters that will be written on the output.
  --   This routine must be invoked before the actual output operation.

  begin
    if n >= right - column
     then new_line(file);
          column := 0;
     else column := column +  n;
    end if;
  end line;

-- AUXILIARIES FOR THE INPUT ROUTINES :

  procedure Build_Number ( file : in file_type;
                           char : in out character; i1,i2 : out integer;
                           ni1,ni2 : out natural; sign : out character ) is
  -- DESCRIPTION :
  --   characters are read from the input and a number is build up;
  --   the result is the number : i1*10^ni2 + i2.
  
  -- ON ENTRY :
  --   file         file type of a file that must be opened for input;
  --   char         first character to be analized.
 
  -- ON RETURN :
  --   char         first character that is not a digit;
  --   i1, i2       digits read;
  --   ni1, ni2     number of digits in i1 and i2;
  --   sign         sign of the number.

    res1,res2 : integer := 0;
    min : boolean := false;
    k1,k2,temp : natural := 0;

  begin
    sign := '+';
    while (char = ' ') loop get(file,char); end loop;
    if (char = '+') or (char = '-')
     then min := (char = '-');
          sign := char;
          get(file,char); 
    end if;
    while (char = ' ') loop get(file,char); end loop;
    loop
      temp := convert(char);
      if temp < 10
       then if k1 < 9
             then res1 := res1*10 + temp;
                  k1 := k1 + 1;
             elsif k2 < 9
                  then res2 := res2*10 + temp;
                       k2 := k2 + 1;
                  else null;  -- skip the rest of the numbers
            end if;
            get(file,char);
       else exit;
      end if;
    end loop;
    if min
     then i1 := -res1; i2 := -res2;
     else i1 := res1;  i2 := res2;
    end if;
    ni1 := k1;
    ni2 := k2;
  end Build_Number;

  procedure Build_Number ( file : in file_type;
                           char : in out character; f : out double_float ) is
  -- DESCRIPTION :
  --  a floating point number is read

    f_int1,f_int2,f_quot1,f_quot2,expo,expo2 : integer := 0;
    f_int,f_quot : double_float := 0.0;
    k1,k2,nq1,nq2,np1,np2,temp : natural := 0;
    sign : character;
    min : boolean;

  begin
    Build_Number(file,char,f_int1,f_int2,np1,np2,sign);
    f_int := double_float(f_int1) * 10.0**np2 + double_float(f_int2);
    min := (sign = '-');
    case char is
      when '.'    => get(file,char);       -- skip the point
                     temp := convert(char);
                     if temp < 10
                      then Build_Number(file,char,f_quot1,f_quot2,nq1,nq2,sign);
                           f_quot := double_float(f_quot1) * 10.0**nq2
                                     + double_float(f_quot2);
                     end if;
                     if char = 'E'
                      then get(file,char); -- skip the 'E'
                           Build_Number(file,char,expo,expo2,k1,k2,sign);
                     end if;
      when 'E'    => if char = 'E'
                      then get(file,char); -- skip the 'E'
                           Build_Number(file,char,expo,expo2,k1,k2,sign);
                     end if;
      when others => null;
    end case; 
    if min
     then if (f_int = 0.0) and (f_quot = 0.0) and (nq1 = 0) and (np1 = 0)
           then f := -1.0;   --  "-x" = -1*x 
           else f := ( f_int - f_quot*10.0**(-nq1-nq2) )*10.0**expo ;
          end if;
     else f := ( f_int + f_quot*10.0**(-nq1-nq2) )*10.0**expo ;
    end if;
  end Build_Number;

  procedure Build_Number ( file : in file_type;
                           char : in out character; c : out double_complex ) is
  -- DESCRIPTION :
  --   a floating point number is read and converted into a complex number;
  --   the number may be the quotient of two floating point numbers
 
    f1,f2 : double_float;

  begin
    Build_Number(file,char,f1);
    if char = '/'
     then get(file,char);            -- skip the '/'
          Build_Number(file,char,f2);
          c := CMPLX(f1/f2);
     else c := CMPLX(f1);
    end if;
  exception
    when numeric_error => raise INFINITE_NUMBER;
  end Build_Number;

  procedure Read_Term ( file : in file_type; char : in out character;
                        n : in natural; termp : in out Poly );
  -- DESCRIPTION :
  --   Reads a term from file, char is the first character of the term.

  procedure Read_Factor ( file : in file_type;
                          char : in out character; n : in natural;
                          d : in out Degrees; pb : in out Poly );
  -- DESCRIPTION :
  --   Reads a factor from file, char is the first character of the factor.

  procedure Read_Factor ( file : in file_type;
                          char : in out character; n : in natural;
                          d : in out Degrees; pb : in out Poly ) is

    sb : symbol;
    i : positive := 1;
    k,ne,ne2 : natural := 0;
    expo,expo2 : integer := 1;
    sign : character;
 
  begin
    sb := (sb'range => ' ');
    while (char = ' ') loop get(file,char); end loop;
    if char = '('
     then get(file,n,pb);
          get(file,char);       -- get a new symbol, skip '('
          return;
    end if;
   -- read the symbol :
    loop
      case char is
        when '+' | '-' | '*' | '^' => exit;
        when delimiter | ' ' | ')' => exit;
        when '('                   => raise ILLEGAL_SYMBOL;
        when others                => sb(i) := char;
                                      i := i+1; get(file,char);
      end case;
    end loop;
   -- check for legality of the symbol :
    if convert(sb(1)) < 10
     then raise ILLEGAL_SYMBOL;
     else for j in 2..3 loop
            case sb(j) is
              when '*' | '+' | '-' | '^' | '/' | ';' | '(' | ')'
                => raise ILLEGAL_SYMBOL;
              when others => null;
            end case;
          end loop;
    end if;
   -- search for the number of the symbol :
    k := Symbol_Table.get(sb);
    if k = 0
     then declare
          begin
            Symbol_Table.add(sb,k);
          exception
            when OVERFLOW_IN_THE_SYMBOL_TABLE => raise OVERFLOW_OF_UNKNOWNS;
          end;
    end if;
    if k > n
     then raise OVERFLOW_OF_UNKNOWNS;
    end if;
   -- read further :
    while (char = ' ') loop get(file,char); end loop;
    if char = '^'
     then get(file,char);                                    -- skip the '^'
          Build_Number(file,char,expo,expo2,ne,ne2,sign);
          d(k) := d(k) + natural(expo);
          while char = ' ' loop get(file,char); end loop;
          if char /= '*'                            -- the case x^2*...
           then return;                             -- end of factor
           else get(file,char);                     -- skip the '*'
          end if; 
     elsif char = '*'
         then get(file,char);
              if char = '*'
               then get(file,char);                 -- the case " x ** expo "
                    Build_Number(file,char,expo,expo2,ne,ne2,sign);
                    d(k) := d(k) + natural(expo);
                    while (char = ' ') loop get(file,char); end loop;
                    if char /= '*'
                     then return;                   -- end of factor
                     else get(file,char);           -- skip the '*'
                    end if;
               else d(k) := d(k) + 1;               -- the case " x * ? "
              end if;
         else -- the case " x ?", with ? /= '*' or ' ' or '^' :
              d(k) := d(k) + 1;
              return;
    end if;
    while (char = ' ') loop get(file,char); end loop;
    if (char = '-') or (char = '+') 
     then return;
    end if;
    if convert(char) < 10
     then -- the case " x * c " or " x ** c * c " :
          Read_Term(file,char,n,pb);
     else -- the case " x * y " :
          Read_Factor(file,char,n,d,pb);
    end if;
  exception
    when ILLEGAL_CHARACTER    => raise ILLEGAL_CHARACTER;
    when ILLEGAL_SYMBOL       => raise ILLEGAL_SYMBOL;
    when ILLEGAL_OPERATION    => raise ILLEGAL_OPERATION;
    when INFINITE_NUMBER      => raise INFINITE_NUMBER;
    when OVERFLOW_OF_UNKNOWNS => raise OVERFLOW_OF_UNKNOWNS;
    when BAD_BRACKET          => raise BAD_BRACKET;
   -- when others               => raise;
  end Read_Factor;
 
  procedure Read_Term ( file : in file_type; char : in out character;
                        n : in natural; termp : in out Poly ) is

    c,c2 : double_complex;
    d : Degrees := new Natural_Vectors.Vector'(1..n => 0);
    pb,res,temp : Poly;
    tmp : Term;

    procedure Collect_Factor_Polynomial is
    begin
      if pb  /= Null_Poly
       then if res = Null_Poly
             then Copy(pb,res); Clear(pb);
             else Mult_Poly(res,pb); Clear(pb);
            end if;
      end if;
    end Collect_Factor_Polynomial;

  begin
    Build_Number(file,char,c);
   
   -- look for 'i' :

    while (char = ' ') loop get(file,char); end loop;
  
    if ( c = CMPLX(0.0) ) and then (char = 'i')
     then -- the case "+ i" :
          c := CMPLX(0.0,1.0); 
          get(file,char);        -- skip 'i'
     elsif ( c = CMPLX(-1.0) ) and then (char = 'i')
         then -- the case "- i" :
              c := CMPLX(0.0,-1.0);
              get(file,char);    -- skip 'i'
         elsif char = '*'
             then -- the case ".. c *.." :
                  while (char = ' ') loop get(file,char); end loop;
                  get(file,char);  -- skip '*'
                  while (char = ' ') loop get(file,char); end loop;
                  if char = 'i'
                   then -- the case ".. c * i.." :
                        c := CMPLX(0.0,REAL_PART(c));
                        get(file,char);    -- skip 'i'
                   else -- the case ".. c * x.." :
                        Read_Factor(file,char,n,d,pb);
                        if pb /= Null_Poly
                         then Clear(res); Copy(pb,res); Clear(pb);
                        end if;
                  end if;
             else -- the case ".. c ?" :
                  -- will be treated in the loop
                  null;
    end if;
 
    loop
      case char is
        when ' '       => get(file,char);
        when '*'       => get(file,char); Read_Factor(file,char,n,d,pb);
                          Collect_Factor_Polynomial;
        when '+' | '-' => if c = CMPLX(0.0)
                           then raise ILLEGAL_CHARACTER;
                           else exit;
                          end if;
        when delimiter => exit;
        when '('       => if c = CMPLX(0.0) or else c = CMPLX(-1.0)
                           then -- the case "+ (" or "- (" :
                                c := CMPLX(0.0);
                                exit;
                           else -- the case "c  (" :
                                raise BAD_BRACKET;
                          end if;
        when ')'       => exit;
        when others    => if c = CMPLX(0.0)
                           then c := CMPLX(1.0);
                                Read_Factor(file,char,n,d,pb);
                           elsif c = CMPLX(-1.0)
                               then Read_Factor(file,char,n,d,pb);
                               else raise ILLEGAL_CHARACTER;
                          end if;
                          Collect_Factor_Polynomial;
      end case;
    end loop;
    tmp.cf := c;
    tmp.dg := d;
    termp := create(tmp);
    if Number_Of_Unknowns(res) > 0
     then Mult_Poly(termp,res); Clear(res);
    end if;
  exception
    when ILLEGAL_CHARACTER    => raise ILLEGAL_CHARACTER;
    when ILLEGAL_SYMBOL       => raise ILLEGAL_SYMBOL;
    when ILLEGAL_OPERATION    => raise ILLEGAL_OPERATION;
    when INFINITE_NUMBER      => raise INFINITE_NUMBER;
    when OVERFLOW_OF_UNKNOWNS => raise OVERFLOW_OF_UNKNOWNS;
    when BAD_BRACKET          => raise BAD_BRACKET;
   -- when others               => raise;
  end Read_Term;

----------------------------------
--    THE INPUT OPERATIONS :    --
----------------------------------

  procedure get ( n : in natural; p : out Poly ) is
  begin
    get(Standard_Input,n,p);
  exception
    when ILLEGAL_CHARACTER    => raise ILLEGAL_CHARACTER;
    when ILLEGAL_SYMBOL       => raise ILLEGAL_SYMBOL;
    when ILLEGAL_OPERATION    => raise ILLEGAL_OPERATION;
    when INFINITE_NUMBER      => raise INFINITE_NUMBER;
    when OVERFLOW_OF_UNKNOWNS => raise OVERFLOW_OF_UNKNOWNS;
    when BAD_BRACKET          => raise BAD_BRACKET;
   -- when others               => raise;
  end get;

  procedure get ( file : in file_type; n : in natural; p : out Poly ) is

    char,oper : character;
    term,res,acc : Poly;

  begin

    if Symbol_Table.Empty
     then Symbol_Table.Init(n);
    end if;

    oper := '+';
    get(file,char);
    while (char = ' ') loop get(file,char); end loop;
    if char = '-'
     then oper := '-';
    end if;
                                    -- the first term can have no sign
    Read_Term(file,char,n,res);     -- therefore read it first
    loop
      case char is
        when ' '       => get(file,char);    -- skip blanks
        when '+' | '-' => oper := char;
                          Read_Term(file,char,n,term);
                          Plus_Poly(res,term); Clear(term);
        when delimiter => exit;
        when '('       => get(file,n,term);
                          case oper is
                            when '+' => Plus_Poly(acc,res); Clear(res);
                                        Copy(term,res);
                            when '-' => Plus_Poly(acc,res);Clear(res);
                                        Copy(term,res); Min_Poly(res);
                            when '*' => Mult_Poly(res,term);
                            when others => raise ILLEGAL_OPERATION;
                          end case;
                          Clear(term);
                          get(file,char);   -- get new character
        when ')'       => exit;
        when '*'       => if res = Null_Poly
                           then raise ILLEGAL_CHARACTER;
                           else -- the case " ) * " :
                                oper := char; get(file,char);  -- skip '*'
                                Read_Term(file,char,n,term);
                                if char /= '('
                                 then case oper is
                                        when '+' => Plus_Poly(res,term);
                                        when '-' => Min_Poly(res,term);
                                        when '*' => Mult_Poly(res,term);
                                        when others => raise ILLEGAL_OPERATION;
                                      end case;
                                end if;
                                Clear(term);
                          end if;
        when others    => raise ILLEGAL_CHARACTER;
      end case;
    end loop;
    p := acc + res;
    Clear(acc); Clear(res);
  exception
    when ILLEGAL_CHARACTER    => raise ILLEGAL_CHARACTER;
    when ILLEGAL_SYMBOL       => raise ILLEGAL_SYMBOL;
    when ILLEGAL_OPERATION    => raise ILLEGAL_OPERATION;
    when INFINITE_NUMBER      => raise INFINITE_NUMBER;
    when OVERFLOW_OF_UNKNOWNS => raise OVERFLOW_OF_UNKNOWNS;
    when BAD_BRACKET          => raise BAD_BRACKET;
   -- when others               => raise;
  end get;
 
  procedure get ( p : out Poly ) is
  begin
    get(Standard_Input,p);
  exception
    when ILLEGAL_CHARACTER    => raise ILLEGAL_CHARACTER;
    when ILLEGAL_SYMBOL       => raise ILLEGAL_SYMBOL;
    when ILLEGAL_OPERATION    => raise ILLEGAL_OPERATION;
    when INFINITE_NUMBER      => raise INFINITE_NUMBER;
    when OVERFLOW_OF_UNKNOWNS => raise OVERFLOW_OF_UNKNOWNS;
    when BAD_BRACKET          => raise BAD_BRACKET;
   -- when others               => raise;
  end get;

  procedure get ( file : in file_type; p : out Poly ) is
    n : natural;
  begin
    get(file,n);
    get(file,n,p);
  exception
    when ILLEGAL_CHARACTER    => raise ILLEGAL_CHARACTER;
    when ILLEGAL_SYMBOL       => raise ILLEGAL_SYMBOL;
    when ILLEGAL_OPERATION    => raise ILLEGAL_OPERATION;
    when INFINITE_NUMBER      => raise INFINITE_NUMBER;
    when OVERFLOW_OF_UNKNOWNS => raise OVERFLOW_OF_UNKNOWNS;
    when BAD_BRACKET          => raise BAD_BRACKET;
   -- when others               => raise;
  end get;

-- AUXILIARIES FOR OUTPUT ROUTINES :

  function Is_Imag ( c : double_complex ) return boolean is
  begin
    return ( REAL_PART(c) = 0.0 );
  end is_imag;

  function Is_Real ( c : double_complex ) return boolean is
  begin
    return ( IMAG_PART(c) = 0.0 );
  end is_real;

  function Is_Integer ( f : double_float ) return boolean is
  begin
    return ( (f - double_float(integer(f))) = 0.0 );
  exception
    when numeric_error => return false;
  end is_integer;
 
  procedure Write_Number ( file : in file_type; i : in integer ) is

  -- DESCRIPTION : 
  --  writes the integer number with only one blank before it

  begin
    for j in 1..8 loop
      if i < integer(10.0**j)
       then line(file,j+1);
            put(file,i,j+1);
            return;
      end if;
    end loop;
    line(file,11); put(file,i);
  end Write_Number;

  procedure Write_Number ( file : in file_type; f : in double_float ) is
  begin
    if is_integer(f)
     then Write_Number(file,integer(f));
     else line(file,21); put(file,f);
    end if;
  end Write_Number;

  procedure Write_Number ( file : in file_type; c : in double_complex ) is
  begin
    if Is_Real(c)
     then Write_Number(file,REAL_PART(c));
     elsif Is_Imag(c)
        then Write_Number(file,IMAG_PART(c));
             line(file,2); put(file,"*i");
        else line(file,1); put(file,"(");
             Write_Number(file,REAL_PART(c));
             if IMAG_PART(c) > 0.0
              then line(file,2); put(file," +");
              else line(file,2); put(file," -");
             end if;
             if IMAG_PART(c) = 1.0
              then line(file,1); put(file,"i");
              elsif IMAG_PART(c) = -1.0
                  then line(file,3); put(file," -i");
                  else Write_Number(file,abs(IMAG_PART(c)));
                       line(file,2); put(file,"*i");
             end if;
             line(file,1); put(file,")");
    end if;
  end Write_Number;

  function Length_Factor ( d,i : natural; standard : boolean;
                           pow : power ) return natural is
  -- DESCRIPTION :
  --   this procedure computes the number of characters needed
  --   for the output of one factor

    l : natural := 0;
    sb : symbol;

  begin
    if standard
     then if i < 10
           then l := l + 2;
           else l := l + 3;
          end if;
     else sb := Symbol_Table.get(i);
          if sb(3) /= ' '
           then l := l + 3;
           elsif sb(2) /= ' '
               then l := l + 2;
               else l := l + 1;
          end if;
    end if;
    if d > 1
     then if pow = '^'
           then l := l + 1;
           else l := l + 2;
          end if;
          if d < 10
           then l := l + 1;
           else l := l + 2;
          end if;
    end if;
    return l;
  end Length_Factor;

  procedure Write_Factor ( file : in file_type; d,i : in natural;
                           standard : in boolean; pow : in power ) is
  -- DESCRIPTION :
  --   Writes the factor corresponding with the ith unknown on file.

    sb : Symbol;

  begin
    if standard
     then put(file,'x');
          if i<10
           then put(file,i,1);
           else put(file,i,2);
          end if;
     else sb := Symbol_Table.get(i); Symbol_Table_io.put(file,sb);
    end if;
    if d > 1
     then if pow = '^'
           then put(file,'^');
           else put(file,"**");
          end if;
          if d < 10
           then put(file,d,1);
           else put(file,d,2);
          end if;
    end if;
  end Write_Factor;

-- THE OUTPUT OPERATIONS :

  procedure put ( n : out natural; p : in Poly; pow : in power := '*' ) is
  begin
    put(Standard_Output,n,p,pow);
  end put;

  procedure put ( file : in file_type; n : out natural; p : in Poly; 
                  pow : in power := '*' ) is

    nn : constant natural := Number_of_Unknowns(p);
    standard : constant boolean := ( Symbol_Table.number < nn );
    first_time : boolean := true;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
 
    -- DESCRIPTION : 
    --   Writes a term is written on file.
 
      passed : boolean;
    begin
      if first_time 
       then first_time := false;
       else if (is_real(t.cf) and then REAL_PART(t.cf) > 0.0)
              or else (is_imag(t.cf) and then IMAG_PART(t.cf) > 0.0)
              or else (not is_real(t.cf) and then not is_imag(t.cf))
             then line(file,1); put(file,'+');
            end if;
      end if;
      if Sum(t.dg) = 0
       then Write_Number(file,t.cf);
       else if ( t.cf - CMPLX(-1.0) ) + CMPLX(1.0) = CMPLX(1.0)
             then line(file,1); put(file,'-');
             elsif ( t.cf - CMPLX(0.0,1.0) ) + CMPLX(1.0) = CMPLX(1.0)
                then line(file,2); put(file,"i*");
                elsif ( t.cf - CMPLX(0.0,-1.0) ) + CMPLX(1.0) = CMPLX(1.0)
                   then line(file,3); put(file,"-i*");
                   elsif (t.cf /= CMPLX(1.0))
                       then Write_Number(file,t.cf); 
                            line(file,1); put(file,'*');
            end if;
            passed := false;
            for i in t.dg'range loop
              if t.dg(i) > 0
               then if passed
                     then line(file,1); put(file,'*');
                     else passed := true;
                    end if;
                    Line(file,Length_Factor(t.dg(i),i,standard,pow));
                    Write_Factor(file,t.dg(i),i,standard,pow);
              end if;
            end loop;
      end if;
      continue := true;
    end Write_Term;
 
    procedure Write_Terms is new Visiting_Iterator (process => Write_Term);

  begin
    init_line;
    n := nn;
    Write_Terms(p);
    line(file,1); put(file,delimiter);
  end put;

  procedure put ( p : in Poly; pow : in power ) is
  begin
    put(Standard_Output,p,pow);
  end put;

  procedure put ( file : in file_type; p : in Poly; pow : in power ) is

    n : natural := Number_of_Unknowns(p);

  begin
    init_line;
    put(file,n,1); put_line(file," ");
    put(file,n,p,pow);
  end put;

  procedure put ( p : in Poly ) is
  begin
    put(Standard_Output,p,'*');
  end put;
 
  procedure put ( file : in file_type; p : in Poly ) is
  begin
    put(file,p,'*');
  end put;

  procedure put_line ( file : in file_type; p : in Poly; pow : in power ) is

    n : constant natural := Number_of_Unknowns(p);
    standard : constant boolean := ( Symbol_Table.Number < n );

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      new_line(file);
      put(file,"+"); Init_Line; Write_Number(file,t.cf);
      if Sum(t.dg) /= 0
       then for i in t.dg'range loop
              if t.dg(i) > 0
               then put(file,'*');
                    Write_Factor(file,t.dg(i),i,standard,pow);
              end if;
            end loop;
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator (process => Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put_line;

  procedure put_line ( p : in Poly; pow : in power ) is
  begin
    put_line(Standard_Output,p,pow);
  end put_line;

  procedure put_line ( p : in Poly ) is
  begin
    put_line(Standard_Output,p,'*');
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly ) is
  begin
    put_line(file,p,'*');
  end put_line;

  procedure Display_Format is

    s : array(1..24) of string(1..65);

  begin
    s( 1):="  A complex multivariate polynomial is denoted as a  sequence  of";
    s( 2):="terms, separated by `+' and terminated by the semicolon `;'.  The";
    s( 3):="brackets '(' and ')' must be used to isolate a sequence of  terms";
    s( 4):="as a factor in a complex multivariate polynomial.                ";
    s( 5):="  A term can be either a coefficient or a  coefficient,  followed";
    s( 6):="by  '*'  and  a  monomial.  If in the latter case the coefficient";
    s( 7):="equals one, then it may be omitted.                              ";
    s( 8):="  A coefficient may be denoted  as  an  integer,  a  rational,  a";
    s( 9):="floating-point or a complex number.                              ";
    s(10):="  A monomial is a sequence of powers of  unknowns,  separated  by";
    s(11):="'*'.   The power operator is represented by '**' or '^'.  It must";
    s(12):="be followed by a positive natural number.  If  the  power  equals";
    s(13):="one, then it may be omitted.                                     ";
    s(14):="  An unknown can be denoted by at most 3 characters.   The  first";
    s(15):="character  must  be a letter and the other two characters must be";
    s(16):="different from '+', '-', '*', '^', '/', ';', '('  and  ')'.   The";
    s(17):="letter i means sqrt(-1), whence it does not represent an unknown.";
    s(18):="The number of unknowns may not  exceed  the  declared  dimension.";
    s(19):="  Some  examples  of  valid  notations  of  complex  multivariate";
    s(20):="polynomials:                                                     ";
    s(21):="  x**2*y + 1/2*z*y**2 - 2*z + y**3 + x - 1E9/-8.E-6* y + 3;      ";
    s(22):="  x^2*y + z*y^2 - 2*z + y^3 + x - y + 3;                         ";
    s(23):="  (1.01 + 2.8*i)*x1**2*x2 + x3**2*x1 - 3*x1 + 2*x2*x3 - 3;       ";
    s(24):="  (x1^2*x2 + x3^2*x1 - 3*x1 + 2*x2*x3 - 3)*x2**2*(x2-1+i);       ";
    for i in s'range loop
      put_line(s(i));
    end loop;
  end Display_Format;

end Complex_Multivariate_Polynomials_io;
