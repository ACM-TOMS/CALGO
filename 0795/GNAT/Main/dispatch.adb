with text_io;                      use text_io;
with Unix_Command_Line;

with mainscal,mainred;
with mainroco,bablroco;
with bablpoco,mainpoco;
with mainsmvc,babldmvc;
with mainphc,bablphc;
with mainvali,bablvali;

procedure Dispatch is

-- BANNERS WITH INFORMATION TO START DIALOGUE WITH USER :

  welcome : constant string :=
    "Welcome to PHC (Polynomial Homotopy Continuation) Version 1.0.";

  author : constant string :=
    "Author is Jan Verschelde (E-mail: na.jverschelde@na-net.ornl.gov).";

  scalban : constant string :=
    "Equation/variable Scaling on polynomial system and solution list.";

  reduban : constant string :=
    "Linear and nonlinear Reduction w.r.t the total degree of the system.";

  rocoban : constant string :=
    "Root counting and Construction of product and polyhedral start systems.";

  mvcban : constant string :=
    "Mixed-Volume Computation by four different lifting strategies.";

  pocoban : constant string :=
    "Polynomial Continuation defined by a homotopy in one parameter.";

  valiban : constant string :=
    "Validation, refinement and purification of computed solution lists.";

-- AVAILABLE OPTIONS :

  options : constant string := "sdpmrvb";
  -- s : scal => scaling of a polynomial system
  -- d : redu => reduction w.r.t. the total degree
  -- p : poco => polynomial continuation
  -- r : roco => root counting methods
  -- m : mvc  => mixed-volume computation
  -- v : vali => validation of solutions
  -- b : batch or black box processing

  option1,option2 : character;
  posi : natural := 0;
  argc : natural := Unix_Command_Line.Number_of_Arguments;

-- UTILITIES FOR PROCESSING THE ARGUMENTS AND OPTIONS :

  function Read_Argument ( k : in natural ) return string is

  -- DESCRIPTION :
  --   Reads the kth argument from the command line.
  --   An argument is a string not proceeded by a `-' character.
  --   The empty string is returned when there is no argument.

    null_string : constant string := "";
    cnt : natural := 0;

  begin
    if argc >= 1
     then for i in 1..argc loop
            declare
              s : constant string := Unix_Command_Line.Argument(i);
            begin
              if s(1) /= '-'
               then cnt := cnt + 1;
                    if k = cnt
                     then return s;
                    end if;
              end if;
            end;
          end loop;
    end if;
    return null_string;
  end Read_Argument;

  function Position ( c : character; s : string ) return natural is

  -- DESCRIPTION :
  --   If the the string contains the character c, then its position
  --   in the string will be returned.  Otherwise s'first-1 will be returned.

  begin
    for i in s'range loop
      if s(i) = c
       then return i;
      end if;
    end loop;
    return s'first-1;
  end Position;

  procedure Read_Next_Option ( pos : in out natural; legal : in string;
                               option : out character ) is

  -- DESCRIPTION :
  --   Reads the next option from the command line arguments.

  -- ON ENTRY :
  --   pos      position in the command line of the last option
  --            that has been read;
  --   legal    string which contains all legal options.

  -- ON RETURN :
  --   pos      the position in the command line of the last option read;
  --   option   is blank when no legal option could be read, otherwise it
  --            contains the next legal option.

    res : character := ' ';
    start : natural := pos+1;

  begin
    if argc >= 1
     then for i in start..argc loop
            declare
              s : constant string := Unix_Command_Line.Argument(i);
            begin
              if s(1) = '-'
               then pos := Position(s(2),legal);
                    if pos >= legal'first
                     then res := legal(pos);
                     else put("The option `"); put(s);
                          put_line("' is not recognised.  Will ignore it...");
                    end if;
              end if;
            end;
            pos := i;
            exit when (res /= ' ');
          end loop;
    end if;
    option := res;
  end Read_Next_Option;

-- DISPATCHING ACCORDING TO OPTIONS :

  procedure Dispatcher ( infile,outfile : in string ) is
  begin
    case option1 is
      when 'b'    => Read_Next_Option(posi,options,option2);
                     case option2 is
                       when 's'    => mainscal(infile,outfile);
                       when 'd'    => mainred(infile,outfile);
                       when 'r'    => bablroco(infile,outfile);
                       when 'm'    => babldmvc(infile,outfile);
                       when 'p'    => bablpoco(infile,outfile);
                       when 'v'    => bablvali(infile,outfile);
                       when others => bablphc(infile,outfile);
                     end case;
      when 's'    => put_line(welcome); put_line(scalban);
                     mainscal(infile,outfile);
      when 'd'    => put_line(welcome); put_line(reduban);
                     mainred(infile,outfile);
      when 'r'    => Read_Next_Option(posi,options,option2);
                     case option2 is
                       when 'b'    => bablroco(infile,outfile);
                       when others => put_line(welcome); put_line(rocoban);
                                      mainroco(infile,outfile);
                     end case;
      when 'm'    => Read_Next_Option(posi,options,option2);
                     case option2 is
                       when 'b'    => babldmvc(infile,outfile);
                       when others => put_line(welcome); put_line(mvcban);
                                      mainsmvc(infile,outfile);
                     end case;
      when 'p'    => Read_Next_Option(posi,options,option2);
                     case option2 is
                       when 'b'    => bablpoco(infile,outfile);
                       when others => put_line(welcome); put_line(pocoban);
                                      mainpoco(infile,outfile);
                     end case;
      when 'v'    => Read_Next_Option(posi,options,option2);
                     case option2 is
                       when 'b'    => bablvali(infile,outfile);
                       when others => put_line(welcome); put_line(valiban);
                                      mainvali(infile,outfile);
                     end case;
      when others => put_line(welcome); mainphc(infile,outfile);
    end case;
  end Dispatcher;
 
begin
  Read_Next_Option(posi,options,option1);
  declare
    nullstring : constant string := "";
    argument : constant string := Read_Argument(1);
    outfile : constant string := Read_Argument(2);
  begin
    if (argument /= "") and then (argument = outfile)
     then new_line; 
          put_line("Input and output file have the same name.");
          put_line("Will ignore output file name...");
          Dispatcher(argument,nullstring);
     else Dispatcher(argument,outfile);
    end if; 
  end;
end Dispatch;
