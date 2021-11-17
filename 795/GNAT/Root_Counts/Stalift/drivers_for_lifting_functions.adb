with integer_io;                        use integer_io;
with Communications_with_User;          use Communications_with_User;
with Numbers_io;                        use Numbers_io;
with Random_Number_Generators;
with Integer_Vectors;
with Integer_Vectors_io;                use Integer_Vectors_io;
with Float_Vectors;
with Float_Vectors_io;                  use Float_Vectors_io;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;

with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Integer_Lifting_Functions;         use Integer_Lifting_Functions;
with Integer_Lifting_Utilities;         use Integer_Lifting_Utilities;
with Float_Lifting_Functions;           use Float_Lifting_Functions;
with Float_Lifting_Utilities;           use Float_Lifting_Utilities;
with Float_Integer_Convertors;          use Float_Integer_Convertors;
with Float_Mixed_Subdivisions_io;

package body Drivers_for_Lifting_Functions is

-- AUXILIARIES :

  procedure Read_Integer_Vector
              ( v : in out Integer_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Give integer for component "); put(i,1); put(" : ");
      Read_Integer(v(i));
    end loop;
  end Read_Integer_Vector;

  procedure Read_Float_Vector
              ( v : in out Float_Vectors.Vector ) is
  begin
    for i in v'range loop
      put("Give float for component "); put(i,1); put(" : ");
      Read_Double_Float(v(i));
    end loop;
  end Read_Float_Vector;

  procedure Set_Integer_Bounds
              ( file : in file_type;
                low,upp : in out Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows the default values for lower and upper bounds and allows
  --   the user to modify these.

    ans : character;
    m : constant natural := low'length;

  begin
    loop
      new_line;
      put_line("Current lower and upper bounds on lifting values");
      put("  lower : "); put(low); new_line;
      put("  upper : "); put(upp); new_line;
      put("Do you want to change these values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      put("Reading "); put(m,1); put(" lower bounds ");
      put_line("for the lifting."); Read_Integer_Vector(low);
      put("Reading "); put(m,1); put(" upper bounds ");
      put_line("for the lifting."); Read_Integer_Vector(upp);
    end loop;
    put(file,"  Lower bounds : "); put(file,low); new_line(file);
    put(file,"  Upper bounds : "); put(file,upp); new_line(file);
  end Set_Integer_Bounds;

  procedure Set_Float_Bounds
              ( file : in file_type;
                low,upp : in out Float_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows the default values for lower and upper bounds and allows
  --   the user to modify these.

    ans : character;
    m : constant natural := low'length;

  begin
    loop
      new_line;
      put_line("Current lower and upper bounds on lifting values");
      put("  lower : "); put(low,2,3,3); new_line;
      put("  upper : "); put(upp,2,3,3); new_line;
      put("Do you want to change these values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      put("Reading "); put(m,1); put(" lower bounds ");
      put_line("for the lifting."); Read_Float_Vector(low);
      put("Reading "); put(m,1); put(" upper bounds ");
      put_line("for the lifting."); Read_Float_Vector(upp);
    end loop;
    put(file,"  Lower bounds : "); put(file,low,2,3,3); new_line(file);
    put(file,"  Upper bounds : "); put(file,upp,2,3,3); new_line(file);
  end Set_Float_Bounds;

  function Read_Integer_Lifting
              ( l : Lists_of_Integer_Vectors.List )
              return Lists_of_Integer_Vectors.List is

    use Lists_of_Integer_Vectors;
    use Integer_Vectors;

    tmp : List := l;
    res,res_last : List;

  begin
    put_line("Give a lifting value for the following points :");
    while not Is_Null(tmp) loop
      declare
        v : constant Vector := Head_Of(tmp).all;
        lv : Vector(v'first..v'last+1);
      begin
        lv(v'range) := v;
        put(v); put(" : "); Read_Integer(lv(lv'last));
        Append(res,res_last,lv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Read_Integer_Lifting;

  function Read_Float_Lifting
              ( l : Lists_of_Float_Vectors.List ) 
              return Lists_of_Float_Vectors.List is

    use Lists_of_Float_Vectors;
    use Float_Vectors;

    tmp : List := l;
    res,res_last : List;

  begin
    put_line("Give a lifting value for the following points :");
    while not Is_Null(tmp) loop
      declare
        v : constant Vector := Head_Of(tmp).all;
        lv : Vector(v'first..v'last+1);
      begin
        lv(v'range) := v;
        put(v,0,0,0); put(" : "); Read_Double_Float(lv(lv'last));
        Append(res,res_last,lv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Read_Float_Lifting;

  procedure Integer_Random_Linear_Lifting
             ( file : in file_type; n : in natural;
               points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lilifu : in out Integer_Vectors_of_Vectors.Link_to_Vector ) is

    low : Integer_Vectors.Vector(points'range) := (points'range => 0);
    upp : Integer_Vectors.Vector(points'range) := Adaptive_Lifting(points);
    ranvec : Integer_Vectors.Vector(1..n);

  begin
    Set_Integer_Bounds(file,low,upp);
    lilifu := new Integer_Vectors_of_Vectors.Vector(points'range);
    for i in lilifu'range loop
      for j in ranvec'range loop
        ranvec(j) := Random_Number_Generators.Random(low(i),upp(i));
      end loop;
      lilifu(i) := new Integer_Vectors.Vector'(ranvec);
    end loop;
    lifted := Linear_Lift(lilifu.all,points);
  end Integer_Random_Linear_Lifting;

  procedure Float_Random_Linear_Lifting
             ( file : in file_type; n : in natural;
               points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
               lilifu : in out Float_Vectors_of_Vectors.Link_to_Vector ) is

    m : constant natural := points'last;
    low : Float_Vectors.Vector(points'range) := (points'range => 0.0);
    upp : Float_Vectors.Vector(points'range) := Adaptive_Lifting(points);
    ranvec : Float_Vectors.Vector(1..n);

  begin
    Set_Float_Bounds(file,low,upp);
    lilifu := new Float_Vectors_of_Vectors.Vector(points'range);
    for i in lilifu'range loop
      ranvec := Random(m,low(i),upp(i));
      lifted(i) := Linear_Lift(points(i),ranvec);
      lilifu(i) := new Float_Vectors.Vector'(ranvec);
    end loop;
  end Float_Random_Linear_Lifting;

  procedure Integer_User_Linear_Lifting
             ( file : in file_type; n : in natural;
               points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lilifu : in out Integer_Vectors_of_Vectors.Link_to_Vector ) is

    lift : Integer_Vectors_of_Vectors.Vector(points'range);

  begin
    for k in lift'range loop
      put("Give "); put(n,1); put(" integers for ");
      put("vector "); put(k,1); put(" : "); get(n,lift(k));
      put(file," vector "); put(file,k,1);
      put(file," : "); put(file,lift(k)); new_line(file);
    end loop;
    lifted := Linear_Lift(lift,points);
    lilifu := new Integer_Vectors_of_Vectors.Vector'(lift);
  end Integer_User_Linear_Lifting;

  procedure Float_User_Linear_Lifting
             ( file : in file_type; n : in natural;
               points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
               lilifu : in out Float_Vectors_of_Vectors.Link_to_Vector ) is

    lift : Float_Vectors_of_Vectors.Vector(points'range);

  begin
    for k in lift'range loop
      put("Give "); put(n,1); put(" floats for ");
      put("vector "); put(k,1); put(" : "); get(n,lift(k));
      put(file," vector "); put(file,k,1);
      put(file," : "); put(file,lift(k)); new_line(file);
      lifted(k) := Linear_Lift(points(k),lift(k).all);
    end loop;
    lilifu := new Float_Vectors_of_Vectors.Vector'(lift);
  end Float_User_Linear_Lifting;

  procedure Integer_Polynomial_Lifting 
            ( file : in file_type;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    lift : Link_to_Poly_Sys;

  begin
    get(lift);
    put(file,lift.all);
    lifted := Polynomial_Lift(lift.all,points);
  end Integer_Polynomial_Lifting;

  procedure Float_Polynomial_Lifting
            ( file : in file_type;
              points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists ) is

    lift : Link_to_Poly_Sys;

  begin
    get(lift);
    put(file,lift.all);
    lifted := Polynomial_Lift(lift.all,points);
  end Float_Polynomial_Lifting;

  procedure Integer_User_Point_Wise_Lifting
            ( file : in file_type;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is
  begin
    for k in points'range loop
      lifted(k) := Read_Integer_Lifting(points(k));
    end loop;
  end Integer_User_Point_Wise_Lifting;

  procedure Float_User_Point_Wise_Lifting
            ( file : in file_type;
              points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists ) is
  begin
    for k in points'range loop
      lifted(k) := Read_Float_Lifting(points(k));
    end loop;
  end Float_User_Point_Wise_Lifting;

  procedure Integer_Random_Point_Wise_Lifting
            ( file : in file_type;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    low : Integer_Vectors.Vector(points'range) := (points'range => 0);
    upp : Integer_Vectors.Vector(points'range) := Adaptive_Lifting(points);

  begin
    Set_Integer_Bounds(file,low,upp);
    lifted := Random_Lift(low,upp,points);
  end Integer_Random_Point_Wise_Lifting;

 procedure Float_Random_Point_Wise_Lifting
            ( file : in file_type;
              points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists ) is

    low : Float_Vectors.Vector(points'range) := (points'range => 0.0);
    upp : Float_Vectors.Vector(points'range) := Adaptive_Lifting(points);

  begin
    Set_Float_Bounds(file,low,upp);
    lifted := Random_Lift(points,low,upp);
  end Float_Random_Point_Wise_Lifting;

  function Menu_for_Lifting_Functions return character is

  -- DESCRIPTION :
  --   Displays the menu and returns the selected lifting function.

    ans : character;

  begin
    put_line("MENU for Lifting Functions :");
    put_line("  Linear     : 1. Random linear vectors will be chosen.");
    put_line("             : 2. You can choose the vectors yourself.");
    put_line("  Polynomial : 3. The system to be solved will be used.");
    put_line("             : 4. You can choose the polynomials yourself.");
    put_line("  Point-wise : 5. For each point a random lifting value.");
    put_line("             : 6. You can choose all lifting values yourself.");
    put("Type a number between 1 and 6 to select lifting : ");
    Ask_Alternative(ans,"123456"); return ans;
  end Menu_for_Lifting_Functions;

  procedure Dispatch_Integer_Lifting
              ( file : in file_type; p : in Poly_Sys; choice : in character;
                points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lilifu : in out Integer_Vectors_of_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Calls the appropriate lifting function to perform the lifting.

    n : constant natural := p'length;

  begin
    new_line(file); put(file,"INTEGER LIFTING FUNCTION : ");
    case choice is
      when '1' => put_line(file,"random linear.");
                  Integer_Random_Linear_Lifting(file,n,points,lifted,lilifu);
      when '2' => put_line(file,"linear provided by user.");
                  Integer_User_Linear_Lifting(file,n,points,lifted,lilifu);
      when '3' => put_line(file,"polynomial system.");
                  lifted := Polynomial_Lift(p,points);
      when '4' => put_line(file,"polynomials provided by user.");
                  Integer_Polynomial_Lifting(file,points,lifted);
      when '5' => put_line(file,"random point-wise.");
                  Integer_Random_Point_Wise_Lifting(file,points,lifted);
      when '6' => put_line(file,"point-wise provided by user.");
                  Integer_User_Point_Wise_Lifting(file,points,lifted);
      when others => null;
    end case;
    new_line(file); put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
    put(file,lifted);
  end Dispatch_Integer_Lifting;

  procedure Dispatch_Float_Lifting
              ( file : in file_type; p : in Poly_Sys; choice : in character;
                points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
                lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
                lilifu : in out Float_Vectors_of_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Calls the appropriate lifting function to perform the lifting.

    n : constant natural := p'length;

  begin
    new_line(file); put(file,"FLOATING-POINT LIFTING FUNCTION : ");
    case choice is
      when '1' => put_line(file,"random linear.");
                  Float_Random_Linear_Lifting(file,n,points,lifted,lilifu);
      when '2' => put_line(file,"linear provided by user.");
                  Float_User_Linear_Lifting(file,n,points,lifted,lilifu);
      when '3' => put_line(file,"polynomial system.");
                  lifted := Polynomial_Lift(p,points);
      when '4' => put_line(file,"polynomials provided by user.");
                  Float_Polynomial_Lifting(file,points,lifted);
      when '5' => put_line(file,"random point-wise.");
                  Float_Random_Point_Wise_Lifting(file,points,lifted);
      when '6' => put_line(file,"point-wise provided by user.");
                  Float_User_Point_Wise_Lifting(file,points,lifted);
      when others => null;
    end case;
    new_line(file); put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
    Float_Mixed_Subdivisions_io.put(file,lifted);
  end Dispatch_Float_Lifting;

-- TARGET ROUTINES :

  procedure Driver_for_Lifting_Functions
            ( file : in file_type; p : in Poly_Sys;
              points : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
              lilifu : in out Integer_Vectors_of_Vectors.Link_to_Vector ) is

    liftfun : character := Menu_for_Lifting_Functions;

  begin
    Dispatch_Integer_Lifting(file,p,liftfun,points,lifted,lilifu);
  end Driver_for_Lifting_Functions;

  procedure Driver_for_Lifting_Functions
            ( file : in file_type; p : in Poly_Sys;
              points : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
              lifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
              lilifu : in out Float_Vectors_of_Vectors.Link_to_Vector ) is

    liftfun : character := Menu_for_Lifting_Functions;

  begin
    Dispatch_Float_Lifting(file,p,liftfun,points,lifted,lilifu);
  end Driver_for_Lifting_Functions;

  procedure Driver_for_Lifting_Functions
              ( file : in file_type; p : in Poly_Sys;
                ipoints : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                fltlif : out boolean;
                fpoints : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
                ilifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                flifted : in out Arrays_of_Float_Vector_Lists.Array_of_Lists;
                ililifu : in out Integer_Vectors_of_Vectors.Link_to_Vector;
                flilifu : in out Float_Vectors_of_Vectors.Link_to_Vector ) is

    liftfun : character := Menu_for_Lifting_Functions;
    ans : character;

  begin
    put("Type i for integer or f for floating-point lifting : ");
    Ask_Alternative(ans,"if"); fltlif := (ans = 'f');
    if ans = 'i'
     then Dispatch_Integer_Lifting(file,p,liftfun,ipoints,ilifted,ililifu);
     else fpoints := Convert(ipoints);
          Dispatch_Float_Lifting(file,p,liftfun,fpoints,flifted,flilifu);
    end if;
  end Driver_for_Lifting_Functions;

end Drivers_for_Lifting_Functions;
