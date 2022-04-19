with integer_io;                        use integer_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with m_Homogeneous_Bezout_Numbers;      use m_Homogeneous_Bezout_Numbers;
with Total_Degree_Start_Systems;        use Total_Degree_Start_Systems;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

with Solutions,Solutions_io;            use Solutions,Solutions_io;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;

with Driver_for_Own_Start_System;
with Drivers_for_m_Homogenization;      use Drivers_for_m_Homogenization;
with Drivers_for_Multi_Homogenization;  use Drivers_for_Multi_Homogenization;
with Drivers_for_Set_Structures;        use Drivers_for_Set_Structures;
with Drivers_for_Symmetric_Set_Structures;
 use Drivers_for_Symmetric_Set_Structures;

with Drivers_for_Implicit_Lifting;      use Drivers_for_Implicit_Lifting;
with Drivers_for_Static_Lifting;        use Drivers_for_Static_Lifting;
with Drivers_for_Dynamic_Lifting;       use Drivers_for_Dynamic_Lifting;
with Drivers_for_Symmetric_Lifting;     use Drivers_for_Symmetric_Lifting;

procedure Driver_for_Root_Counts 
             ( file : in file_type; p,q : in out Poly_Sys;
               own : in boolean;
               qsols : in out Solution_List; roco : out natural ) is

  timer : timing_widget;
  rc : natural := Total_Degree(p);
  lpos : List;
  choice : string(1..2) := "  ";
  method,ans : character := 'y';
  noqsols : natural := 0;

  procedure High_Total_Degree is
  begin
    for i in p'range loop
      put(Degree(p(i)),1);    put(file,Degree(p(i)),1);
      exit when i = p'last;
      put("*");               put(file,"*");
    end loop;
    new_line;
    put_line("  this is higher than my largest integer.  Be careful...");
  end High_Total_Degree;

  procedure Display_Menu ( rc : in natural ) is

    m : array(0..9) of string(1..66);

  begin
    new_line;
    put_line("MENU with ROOT COUNTS and Methods to Construct START SYSTEMS :");
    put("  0. exit - current start system is ");
    if Is_Null(qsols)
     then put("based on total degree : "); put(rc,1); new_line;
     else case method is
            when '1' => put("based on m-homogenization : ");
            when '2' => put("based on multi-homogenization : ");
            when '3' => put("based on set structure : ");
            when '4' => put("based on symmetric set structure : ");
            when '5' => put("based on Bezout and BKK Bound : ");
            when '6' => put("based on static mixed-volume computation : ");
            when '7' => put("based on dynamic mixed-volume computation : ");
            when '8' => put("based on symmetric mixed-volume computation : ");
            when '9' => put("your start system : ");
            when others => put("based on total degree");
          end case;
          put(rc,1); new_line;
    end if;
    m(0):="PRODUCT HOMOTOPIES based on DEGREES ------------------------------";
    m(1):="  1. m-homogeneous Bezout number                   (one partition)";
    m(2):="  2. multi-homogeneous Bezout number             (many partitions)";
    m(3):="  3. generalized Bezout number                     (set structure)";
    m(4):="  4. symmetric generalized Bezout number (symmetric set structure)";
    m(5):="POLYHEDRAL HOMOTOPIES based on NEWTON POLYTOPES ------------------";
    m(6):="  5. combination between Bezout and BKK Bound   (implicit lifting)";
    m(7):="  6. mixed-volume computation                     (static lifting)";
    m(8):="  7. incremental mixed-volume computation        (dynamic lifting)";
    m(9):="  8. symmetric mixed-volume computation        (symmetric lifting)";
    for i in m'range loop
      put_line(m(i));
    end loop;
    if own
     then put_line
         ("START SYSTEM DEFINED BY USER -------------------------------------");
          put_line("  9. you can give your own start system");
     else put_line
         ("------------------------------------------------------------------");
    end if;
  end Display_Menu;

  procedure Display_Info ( method : character ) is

  -- DESCRIPTION :
  --   Displays the information that corresponds with the current method.

  begin
    new_line;
    case method is
      when '0' => Total_Degree_Info;
      when '1' => m_Homogenization_Info;
      when '2' => Multi_Homogenization_Info;
      when '3' => Set_Structure_Info;
      when '4' => Symmetric_Set_Structure_Info;
      when '5' => Implicit_Lifting_Info;
      when '6' => Static_Lifting_Info;
      when '7' => Dynamic_Lifting_Info;
      when '8' => Symmetric_Lifting_Info;
      when others => put_line("No information available.");
    end case;
    new_line;
  end Display_Info;

  procedure Apply_Method ( method : character ) is

  -- DESCRIPTION :
  --   Applies the root count that corresponds with the current method.

  begin
    case method is
      when '0' => Start_System(p,q,qsols);
      when '1' => Driver_for_m_Homogenization(file,p,rc,q,qsols);
      when '2' => Driver_for_Multi_Homogenization(file,p,rc,q,qsols);
      when '3' => Driver_for_Set_Structure(file,p,rc,lpos,q,qsols);
      when '4' => Driver_for_Symmetric_Random_Product_Systems
                                             (file,p,q,qsols,rc,lpos);
      when '5' => Driver_for_Mixture_Bezout_BKK(file,p,false,q,qsols,rc);
      when '6' => Driver_for_Mixed_Volume_Computation(file,p,false,q,qsols,rc);
      when '7' => Driver_for_Dynamic_Mixed_Volume_Computation
                                                     (file,p,false,q,qsols,rc);
      when '8' => Driver_for_Symmetric_Mixed_Volume_Computation
                                                     (file,p,false,q,qsols,rc);
      when '9' => Driver_for_Own_Start_System(file,p,q,qsols);
      when others => null;
    end case;
  end Apply_Method;

begin
  new_line(file); put_line(file,"ROOT COUNTS :"); new_line(file);
  put(file,"total degree : ");
  if rc > 0
   then put(file,rc,1); -- put(rc,1);
   else High_Total_Degree;
  end if;
  new_line(file);
  loop
    Display_Menu(rc);
    if own
     then put("Type a number between 0 and 9," 
                & " eventually preceded by i for info : ");
          Ask_Alternative(choice,"0123456789",'i');
     else put("Type a number between 0 and 8," 
                & " eventually preceded by i for info : ");
          Ask_Alternative(choice,"012345678",'i');
    end if;
    if choice(1) = 'i'
     then method := choice(2);
          Display_Info(method);
          put("Do you want to apply this root count ? (y/n) ");
          Ask_Yes_or_No(ans);
     else method := choice(1);
    end if;
    if ans = 'y'
     then Apply_Method(method);
          noqsols := Length_Of(qsols);
          if method /= '0'
           then new_line;
                put("The current root count equals "); put(rc,1); put_line(".");
                if noqsols /= 0
                 then put("The number of start solutions equals ");
                      put(noqsols,1); put_line(".");
                end if;
                put("Do you want to perform more root counting ? (y/n) ");
                Ask_Yes_or_No(ans);
           else ans := 'n';
          end if;
     else ans := 'y';
    end if;
    exit when ans /= 'y';
  end loop;
  roco := rc;
  Clear(lpos);
end Driver_for_Root_Counts;
