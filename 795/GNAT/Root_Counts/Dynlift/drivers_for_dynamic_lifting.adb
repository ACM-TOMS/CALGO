with text_io,integer_io,Numbers_io;      use text_io,integer_io,Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;

with Integer_Vectors,Float_Vectors;      use Integer_Vectors;
with Integer_Vectors_io;                 use Integer_Vectors_io;
with Integer_Vectors_of_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Lists_of_Vectors_Utilities;         use Lists_of_Vectors_Utilities;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;

with Drivers_for_Vertex_Points;          use Drivers_for_Vertex_Points;
with Simplices,Triangulations;           use Simplices,Triangulations;
with Triangulations_io;                  use Triangulations_io;
with Dynamic_Triangulations;             use Dynamic_Triangulations;
with Cayley_Trick;                       use Cayley_Trick;
with Driver_for_Minkowski_Polynomials;

with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;
with Triangulations_and_Subdivisions;    use Triangulations_and_Subdivisions;
with Dynamic_Mixed_Subdivisions;         use Dynamic_Mixed_Subdivisions;
with Dynamic_Polyhedral_Continuation;    use Dynamic_Polyhedral_Continuation;
with Driver_for_Polyhedral_Continuation;

with Power_Lists;                        use Power_Lists;
with Complex_Polynomial_Systems_io;      use Complex_Polynomial_Systems_io;
with Solutions_io;                       use Solutions_io;

with Pruning_Statistics;

package body Drivers_for_Dynamic_Lifting is

  procedure Dynamic_Lifting_Info is

    i : array(1..6) of string(1..65);

  begin
    i(1):="  Dynamic  lifting  can  be  used  to   compute   mixed   volumes";
    i(2):="incrementally,  i.e.:  by  adding  the  points  repeatedly to the";
    i(3):="already constructed subdivision.  This method  works  efficiently";
    i(4):="when  all  Newton polytopes are (almost) equal.  The Cayley trick";
    i(5):="is implemented by means of dynamic lifting.  This trick  computes";
    i(6):="all cells in a mixed subdivision.                                ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Dynamic_Lifting_Info;

  procedure Driver_for_Dynamic_Mixed_Volume_Computation 
                ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                  q : out Poly_Sys; qsols : out Solution_List;
                  mv : out natural ) is

    welcome : constant string := "Mixed-Volume Computation by Dynamic Lifting";

  -- GLOBAL VARIABLES :

    supports : Array_of_Lists(p'range);
    n : natural := p'last;
    timer : timing_widget;
    r,max : natural;
    mix,perms : Link_to_Vector;
    ans : character;
    permp,qq : Poly_Sys(p'range);
    qqsols : Solution_List;
    subfile,solsft,qft : file_type;
    vol : natural := 0;
    mixsub : Mixed_Subdivision;

  -- GLOBAL SWITCHES :

    verpts     : boolean;  -- if the set of vertex points is computed
    order      : boolean;  -- process points in fixed instead of random order
    inter      : boolean;  -- if interior points are possible
    conmv      : boolean;  -- if checks on zero contributions have to be made
    caytrick   : boolean;  -- if the Cayley trick has to be applied
    reportnew  : boolean;  -- if the new cells have to be reported 
    reportflat : boolean;  -- if before flattening, reporting has to be done
    subonfile  : boolean;  -- put the subdivision on separate file
    tosolve    : boolean;  -- if the system needs to be solved
    contrep    : boolean;  -- if intermediate output during continuation
    ranstart   : boolean;  -- if random coefficient start system

    minkpoly   : natural;  -- 0 : no; 1 : only poly, > 1 : all subdivisions

  -- EXTRACT THE ADDITIONAL POINTS :

    function Is_In_Lifted ( pt : Link_to_Vector; lifted : List )
                          return boolean is

    -- DESCRIPTION :
    --   Returns true if the point is in the lifted list.

      tmp : List := lifted;
      lpt : Link_to_Vector;

    begin
      while not Is_Null(tmp) loop
        lpt := Head_Of(tmp);
        if pt(pt'range) = lpt(pt'range)
         then return true;
         else tmp := Tail_Of(tmp);
        end if;
      end loop;
      return false;
    end Is_In_Lifted;

    function Difference ( supp,liftsupp : in List ) return List is

      res,res_last : List;
      tmp : List := supp;
      pt : Link_to_Vector;
 
    begin
      tmp := supp;
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        if not Is_In_Lifted(pt,liftsupp)
         then Append(res,res_last,pt.all);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      return res;
    end Difference;

    function Difference ( supp,liftsupp : in Array_of_Lists )
                        return Array_of_Lists is

    -- DESCRIPTION :
    --   Returns a tuple of point lists, made of points in supp
    --   that do not belong to the corresponding lifted supports.

      res : Array_of_Lists(supp'range);

    begin
      for i in supp'range loop
        res(i) := Difference(supp(i),liftsupp(i));
      end loop;
      return res;
    end Difference;

  -- DETERMINING THE ORDER OF PROCESSING THE POINTS :

    function Determine_Order ( l : List ) return List is

    -- DESCRIPTION :
    --   Interactive ordering of the points in the list.
    --   This function displays all points and asks the user for a position.

      len : constant natural := Length_Of(l);
      pos : vector(1..len);
      res : List;

      function Read_New_Positions
                  ( l : List; length : natural ) return vector is

      -- DESCRIPTION :
      --   Lists all points in the lists and prompts for a new position.
      --   Returns the position vector.

        newpos : vector(1..length);
        tmp : List := l;
        cnt : natural := 0;

      begin
        put("There are "); put(length,1); put_line(" points to order.");
        put_line("Give for each separate point its new position :");
        while not Is_Null(tmp) loop
          cnt := cnt + 1;
          loop
            put(Head_Of(tmp)); put(" : ");
            Read_Natural(newpos(cnt));
            exit when (newpos(cnt) >= 1) and (newpos(cnt) <= length);
            put("New position out of range 1.."); put(length,1);
            put_line(".  Please try again.");
          end loop;
          tmp := Tail_Of(tmp);
        end loop;
        return newpos;
      end Read_New_Positions;

      function Get ( l : List; pos : natural ) return Link_to_Vector is

      -- DESCRIPTION :
      --   Returns the point on the indicated position in the list l.

        tmp : List := l;
        res : Link_to_Vector;

      begin
        if not Is_Null(l)
         then for i in 1..(pos-1) loop
                tmp := Tail_Of(tmp);
                exit when Is_Null(tmp);
              end loop;
              if not Is_Null(tmp)
               then res := Head_Of(tmp);
              end if;
        end if;
        return res;
      end Get;

      function Sort ( l : in List; pos : in vector ) return List is

      -- DESCRIPTION :
      --   Sorts the given list according to the given position vector:
      --   pos(i) determines the new position of the ith point in the list.
      --   If the returning list is empty, then the position vector was
      --   not a permutation.

        empty,res,res_last : List;
        index : natural;

      begin
        for i in pos'range loop   -- search index : pos(index) = i
        index := 0;
          for j in pos'range loop
            if pos(j) = i
             then index := j;
            end if;
            exit when (index /= 0);
          end loop;
          exit when (index = 0);
          Append(res,res_last,get(l,index));  -- append the vector
        end loop;
        if index = 0
         then return empty;
         else return res;
        end if;
      end Sort;

    begin
      if Is_Null(l)
       then return l;
       else loop
              pos := Read_New_Positions(l,len);
              res := Sort(l,pos);
              exit when not Is_Null(res);
              put_line("The given position vector was not a permutation.");
              put_line("Please try again...");
            end loop;
            return res;
      end if;
    end Determine_Order;

    procedure Determine_Processing_Order
                 ( supports : in out Array_of_Lists; mix : in Link_to_Vector;
                   fixed : out boolean ) is

      choice : character;
      cnt : natural;

    begin
      new_line;
      put_line("MENU for the Order of the points to add : ");
      put_line("  1. fixed order, given by the monomial ordering");
      put_line("  2. random order, generated by the algorithm");
      put_line("  3. interactively defined by you");
      put("Type 1,2, or 3 : "); Ask_Alternative(choice,"123");
      case choice is
        when '1'    => fixed := true;
        when '2'    => fixed := false;
        when others => fixed := true;
                       cnt := supports'first;
                       for i in mix'range loop
                         supports(cnt) := Determine_Order(supports(cnt));
                         cnt := cnt + mix(i);
                       end loop;
      end case;
    end Determine_Processing_Order;

  -- INSTANTIATIONS OF THE GENERICS :

    procedure Report_New_Simplices
                     ( t : in Triangulation; point : in Vector ) is

    -- DESCRIPTION :
    --   Writes the new simplices on file and computes their volume.

      v : natural;

    begin
      new_line(file);
      put(file,"The new simplices by adding "); put(file,point);
      put_line(file," : ");
      put(file,n,t,v);
      put(file," with volume addition : ");
      put(file,vol,1); put(file," + "); put(file,v,1);
      vol := vol + v; put(file," = "); put(file,vol,1); put_line(file,".");
    end Report_New_Simplices;
    procedure R_Dynamic_Lifting is
      new Dynamic_Triangulations.Dynamic_Lifting_with_New(Report_New_Simplices);

    procedure Collect_Flattening ( t : in Triangulation; l : List ) is

    -- DESCRIPTION :
    --   Updates the subdivision mixsub with the flattened cells.
    --   The triangulation on entry contains the whole triangulation,
    --   not just the new cells.

      cells : Mixed_Subdivision;

    begin
      if Is_Null(mixsub)
       then cells := Deep_Create(n,t);
       else cells := Non_Flat_Deep_Create(n,t);
            Construct(Head_Of(mixsub),cells);
      end if;
      Flatten(cells);
      mixsub := cells;
    end Collect_Flattening;

    procedure Report_Flattening
                    ( t : in Triangulation; l : in List ) is

    -- DESCRIPTION :
    --   Writes the list of lifted points and the triangulation on file
    --   and updates the mixed subdivision.

    begin
      new_line(file);
      put_line(file,"The list of lifted points before flattening : ");
      put(file,l);
      new_line(file);
      put_line(file,"The triangulation before flattening : ");
      put(file,n,t,vol);
      put(file," with volume "); put(file,vol,1); put_line(file,".");
      Collect_Flattening(t,l);
    end Report_Flattening;
    procedure C_Dynamic_Lifting is
      new Dynamic_Triangulations.Dynamic_Lifting_with_Flat(Collect_Flattening);
    procedure F_Dynamic_Lifting is
      new Dynamic_Triangulations.Dynamic_Lifting_with_Flat(Report_Flattening);
    procedure FR_Dynamic_Lifting is
      new Dynamic_Triangulations.Dynamic_Lifting_with_Flat_and_New
           ( Before_Flattening     => Report_Flattening,
             Process_New_Simplices => Report_New_Simplices);

    procedure Report_New_Cells 
                    ( mixsub : in out Mixed_Subdivision;
                      i : in natural; point : in Vector ) is

    -- DESCRIPTION :
    --   Writes the new mixed cells on file and computes the mixed volume.

      v : natural;

    begin
      if not Is_Null(mixsub)
       then
         new_line(file);
         put(file,"The new mixed cells by adding "); put(file,point);
         new_line(file);
         put(file," to the "); put(file,i,1); put_line(file,"th component : ");
         put(file,n,mix.all,mixsub,v);
         put(file," with volume addition : ");
         put(file,vol,1); put(file," + "); put(file,v,1);
         vol := vol + v; put(file," = "); put(file,vol,1); new_line(file);
      end if;
    end Report_New_Cells;
    procedure R_Dynamic_Cayley is
      new Cayley_Trick.Dynamic_Cayley_with_New(Report_New_Cells);
    procedure Rt_Dynamic_Cayley is
      new Cayley_Trick.Dynamic_Cayley_with_Newt(Report_New_Cells);
    procedure R_Dynamic_Lifting is
      new Dynamic_Mixed_Subdivisions.Dynamic_Lifting_with_New(Report_New_Cells);

    procedure Report_Flattening
                 ( mixsub : in out Mixed_Subdivision;
                   lifted : in Array_of_Lists ) is

    -- DESCRIPTION :
    --   Writes the list of lifted points and the subdivision on file.

    begin
      new_line(file);
      put_line(file,"The list of lifted points before flattening : ");
      for i in lifted'range loop
        put(file," points of "); put(file,i,1);
        put_line(file,"th component : ");
        put(file,lifted(i));
      end loop;
      new_line(file);
      put_line(file,"The mixed subdivision before flattening : ");
      put(file,n,mix.all,mixsub,vol);
      put(file," with volume "); put(file,vol,1); put_line(file,".");
    end Report_Flattening;
    procedure F_Dynamic_Cayley is
      new Cayley_Trick.Dynamic_Cayley_with_Flat(Report_Flattening);
    procedure Ft_Dynamic_Cayley is
      new Cayley_Trick.Dynamic_Cayley_with_Flatt(Report_Flattening);
    procedure FR_Dynamic_Cayley is
      new Cayley_Trick.Dynamic_Cayley_with_Flat_and_New
            (Before_Flattening => Report_Flattening,
             Process_New_Cells => Report_New_Cells);
    procedure FRt_Dynamic_Cayley is
      new Cayley_Trick.Dynamic_Cayley_with_Flat_and_Newt
            (Before_Flattening => Report_Flattening,
             Process_New_Cells => Report_New_Cells);

    procedure Report_Flattening
                 ( mixsub : in out Mixed_Subdivision;
                   fs : in Face_Structures ) is

    -- DESCRIPTION :
    --   Writes the list of lifted points and the subdivision on file.

    begin
      new_line(file);
      put_line(file,"The lists of lifted points before flattening : ");
      for i in fs'range loop
        put(file," points of "); put(file,i,1);
        put_line(file,"th component : ");
        put(file,fs(i).l);
      end loop;
      new_line(file);
      put_line(file,"The mixed subdivision before flattening : ");
      put(file,n,mix.all,mixsub,vol);
      put(file," with volume "); put(file,vol,1); put_line(file,".");
    end Report_Flattening;
    procedure F_Dynamic_Lifting is
      new Dynamic_Mixed_Subdivisions.Dynamic_Lifting_with_Flat
            (Report_Flattening);
    procedure FR_Dynamic_Lifting is
      new Dynamic_Mixed_Subdivisions.Dynamic_Lifting_with_Flat_and_New
            (Before_Flattening => Report_Flattening,
             Process_New_Cells => Report_New_Cells);

  -- MAIN CONSTRUCTORS :

    procedure Compute_Triangulation is

    -- DESCRIPTION :
    --   Application of the dynamic lifting algorithm 
    --   to compute a triangulation of one polytope.

      t : Triangulation;
      support,lifted,lifted_last : List;
      arlifted : Array_of_Lists(mix'range);

    begin
      support := supports(supports'first);
      if verpts
       then Vertex_Points(file,support);
      end if;
      new_line(file);
      put_line(file,"CREATION OF THE TRIANGULATION :");
      new_line(file);
      tstart(timer);
      if reportnew
       then 
         if reportflat
          then FR_Dynamic_Lifting(support,order,inter,max,lifted,lifted_last,t);
          else R_Dynamic_Lifting(support,order,inter,max,lifted,lifted_last,t);
         end if;
       elsif reportflat
           then 
             F_Dynamic_Lifting(support,order,inter,max,lifted,lifted_last,t);
           elsif subonfile
               then
                 C_Dynamic_Lifting
                                (support,order,inter,max,lifted,lifted_last,t);
               else
                 Dynamic_Lifting(support,order,inter,max,lifted,lifted_last,t);
      end if;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"computing the triangulation");
      new_line(file);
      put_line(file,"THE LIFTED SUPPORTS :"); new_line(file);
      put(file,lifted);
      new_line(file);
      put_line(file,"THE TRIANGULATION :"); new_line(file);
      tstart(timer);
      put(file,n,t,vol);
      tstop(timer);
      new_line(file);
      put(file,"The volume : "); put(file,vol,1); new_line(file);
      new_line(file);
      print_times(file,timer,"computing the volume");
      if subonfile
       then if Is_Null(mixsub)
             then put(subfile,n,1); new_line(subfile);
                  put(subfile,1,1); new_line(subfile); -- type of mixture
                  put(subfile,n,t);
             else declare
                    lastcells : Mixed_Subdivision := Non_Flat_Deep_Create(n,t);
                  begin
                    Construct(Head_Of(mixsub),lastcells);
                    mixsub := lastcells;
                    put(subfile,n,mix.all,mixsub);
                  end;
            end if;
            Close(subfile);
      end if;
      mv := vol;
    end Compute_Triangulation;

    procedure Compute_Cayley_Triangulation is

    -- DESCRIPTION :
    --   Application of the dynamic lifting algorithm to compute a mixed 
    --   subdivision of a tuple of polytopes by means of the Cayley trick.

      supp,lifted : Array_of_Lists(1..r);
      t : Triangulation;
      numtri : natural;
  
    begin
      if verpts
       then Vertex_Points(file,mix,supports);
      end if;
      supp := Typed_Lists(mix.all,supports);
      new_line(file);
      put_line(file,"CREATION OF THE MIXED SUBDIVISION :");
      new_line(file);
      tstart(timer);
      if reportnew
       then 
         if reportflat
          then 
            if minkpoly > 0
             then FRt_Dynamic_Cayley(n,mix.all,supp,order,inter,max,lifted,t);
             else FR_Dynamic_Cayley(n,mix.all,supp,order,inter,max,
                                    lifted,mixsub,numtri);
            end if;
          else
            if minkpoly > 0
             then Rt_Dynamic_Cayley(n,mix.all,supp,order,inter,max,lifted,t);
             else R_Dynamic_Cayley(n,mix.all,supp,order,inter,max,lifted,
                                   mixsub,numtri);
            end if;
         end if;
       elsif reportflat
           then
             if minkpoly > 0
              then Ft_Dynamic_Cayley(n,mix.all,supp,order,inter,max,lifted,t);
              else F_Dynamic_Cayley(n,mix.all,supp,order,inter,max,
                                    lifted,mixsub,numtri);
             end if;
           else
             if minkpoly > 0
              then Dynamic_Cayley(n,mix.all,supp,order,inter,max,lifted,t);
              else Dynamic_Cayley(n,mix.all,supp,order,inter,max,
                                  lifted,mixsub,numtri);
             end if;
      end if;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Computing the mixed subdivision");
      new_line(file);
      put_line(file,"THE LIFTED SUPPORTS :");
      new_line(file);
      put(file,lifted);
      if minkpoly > 0
       then declare
              alltri : boolean := (minkpoly > 1);
            begin
              Driver_for_Minkowski_Polynomials(file,n,mix.all,t,alltri,mixsub);
              numtri := Length_Of(t);
            end;
      end if;
      new_line(file);
      put_line(file,"THE MIXED SUBDIVISION :");
      new_line(file);
      tstart(timer);
      put(file,n,mix.all,mixsub,vol);
      tstop(timer);
      new_line(file);
      put(file,"The mixed volume equals : "); put(file,vol,1);
      new_line(file);
      put(file,"Number of cells in auxiliary triangulation : ");
      put(file,numtri,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Computing the mixed volume");
      if subonfile
       then put(subfile,n,mix.all,mixsub);
            Close(subfile);
      end if;
      mv := vol;
    end Compute_Cayley_Triangulation;

    procedure Report_Results
                ( file : in file_type; n : in natural; mix : in Link_to_Vector;
                  mixsub : in out Mixed_Subdivision;
                  fs : in Face_Structures ) is
    begin
      new_line(file);
      put_line(file,"THE LIFTED SUPPORTS :");
      new_line(file);
      for i in fs'range loop
        put(file,fs(i).l);  new_line(file);
      end loop;
      put_line(file,"THE MIXED SUBDIVISION :");
      new_line(file);
      tstart(timer);
      if r = 1
       then put(file,n,fs(fs'first).t,vol);
       else put(file,n,mix.all,mixsub,vol);
      end if;
      tstop(timer);
      new_line(file);
      put(file,"The mixed volume equals : "); put(file,vol,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Computing the mixed volume");
      if subonfile
       then put(subfile,n,mix.all,mixsub);
            Close(subfile);
      end if;
      mv := vol;
    end Report_Results;
  
    procedure Compute_Mixed_Subdivision is

    -- DESCRIPTION :
    --   Application of the dynamic lifting algorithm
    --   to compute a mixed subdivision of a tuple of polytopes.
  
      supp,lifted : Array_of_Lists(1..r);
      fs : Face_Structures(1..r);
      nbsucc,nbfail : Float_Vectors.Vector(1..r) := (1..r => 0.0);

    begin
      if verpts
       then Vertex_Points(file,mix,supports);
      end if;
      supp := Typed_Lists(mix.all,supports);
      new_line(file);
      put_line(file,"CREATION OF THE MIXED SUBDIVISION :");
      new_line(file);
      tstart(timer);
      if reportnew
       then
         if reportflat
          then FR_Dynamic_Lifting
                 (n,mix.all,supp,order,inter,conmv,max,mixsub,fs,nbsucc,nbfail);
          else R_Dynamic_Lifting
                 (n,mix.all,supp,order,inter,conmv,max,mixsub,fs,nbsucc,nbfail);
         end if;
       elsif reportflat
           then F_Dynamic_Lifting
                 (n,mix.all,supp,order,inter,conmv,max,mixsub,fs,nbsucc,nbfail);
           else Dynamic_Lifting 
                 (n,mix.all,supp,order,inter,conmv,max,mixsub,fs,nbsucc,nbfail);
      end if;
      tstop(timer);
      Pruning_Statistics(file,nbsucc,nbfail);
      new_line(file);
      print_times(file,timer,"Computing the mixed subdivision");
      Report_Results(file,n,mix,mixsub,fs);
    end Compute_Mixed_Subdivision;

    procedure Solve_Coefficient_System is

    -- DESCRIPTION :
    --   Application of the dynamic lifting algorithm
    --   to compute a mixed subdivision of a tuple of polytopes and
    --   to solve a start system, with randomized coefficients.

      supp : Array_of_Lists(1..r);
      fs : Face_Structures(1..r);
      lifted : Array_of_Lists(1..r);
      numtri : natural := 0;
      lif,lif_last : List;
      nbsucc,nbfail : Float_Vectors.Vector(1..r) := (1..r => 0.0);

    begin
      if verpts
       then Vertex_Points(file,mix,supports);
            qq := Select_Terms(qq,supports);
      end if;
      supp := Typed_Lists(mix.all,supports);
      new_line(file);
      put_line(file,"SOLVING THE RANDOM COEFFICIENT SYSTEM :");
      new_line(file);
      tstart(timer);
      if mix'last = mix'first
       then
         Dynamic_Unmixed_Solve
               (file,n,supp(supp'first),order,inter,max,fs(fs'first).l,
                fs(fs'first).last,fs(fs'first).t,qq,qqsols);
       else
         if caytrick
          then
            Dynamic_Cayley_Solve(file,n,mix.all,supp,order,inter,max,lifted,
                                 mixsub,numtri,qq,qqsols);
            for i in lifted'range loop
              fs(i).l := lifted(i);
            end loop;
          else
            Dynamic_Mixed_Solve
                     (file,n,mix.all,supp,order,inter,conmv,max,mixsub,fs,
                      nbsucc,nbfail,qq,qqsols);
         end if;
      end if;
      tstop(timer);
      if mix'last > mix'first and not caytrick
       then Pruning_Statistics(file,nbsucc,nbfail);
      end if;
      new_line(file);
      print_times(file,timer,"Computing the solution list");
      Report_Results(file,n,mix,mixsub,fs);
      q := qq; qsols := qqsols;
      if not ranstart
       then put(solsft,qqsols);
            Close(solsft);
      end if;
      if ranstart
       then new_line(qft); put_line(qft,"THE SOLUTIONS :"); new_line(qft);
            put(qft,Length_Of(qqsols),n,qqsols);
            Close(qft);
      end if;
    end Solve_Coefficient_System;

  begin
    new_line; put_line(welcome);
   -- READING GENERAL INPUT INFORMATION :
    supports := Construct_Power_Lists(p);
    new_line;
    put("Do you want to enforce a type of mixture ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y'
     then Compute_Mixture(supports,mix,perms); r := mix'last;
     else put("Give number of different supports : "); Read_Natural(r);
          put("Give vector of occurrences : "); get(r,mix);
          perms := new Vector(1..n);
          for i in perms'range loop
            perms(i) := i;
          end loop;
    end if;
   -- DETERMINE THE GLOBAL SWITCHES :
    put("Do you first want to extract the vertex points ? (y/n) ");
    Ask_Yes_or_No(ans);
    verpts := (ans = 'y');
    inter := not verpts;
    put("Do you have a maximum lifting value ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put("  Give the maximum lifting value : ");
          Read_Positive(max);
     else max := 0;
    end if;
    Determine_Processing_Order(supports,mix,order);
    if (r > 1) 
     then new_line;
          put_line("MENU for Cayley trick : ");
          put_line("  0. No Cayley trick, pruning for mixed cells.");
          put_line("  1. Cayley trick : auxiliary triangulation.");
          put_line("  2. Cayley trick with Minkowski-polynomial.");
          put_line("  3. Cayley trick with all subdivisions.");
          put("Type 0,1,2, or 3 : ");
          Ask_Alternative(ans,"0123");
          caytrick := not (ans = '0');
          case ans is
            when '2' => minkpoly := 1;
            when '3' => minkpoly := 2;
            when others => minkpoly := 0;
          end case;
          if not caytrick
           then put("Do you want online checks on zero contributions ? (y/n) ");
                Ask_Yes_or_No(ans);
                conmv := (ans = 'y');
           else conmv := false;
          end if;
     else caytrick := false; conmv := false;
    end if;
    put("Do you want to have the subdivision on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then subonfile := true;
          put_line("Reading the name of the file.");
          Read_Name_and_Create_File(subfile);
     else subonfile := false;
    end if;
    new_line;
    put("Are the cells to be written on file, during computation ? (y/n) ");
    Ask_Yes_or_No(ans);
    reportnew := (ans = 'y');
    put("Are the cells to be written on file, before flattening ? (y/n) ");
    Ask_Yes_or_No(ans);
    reportflat := (ans = 'y');
    permp := Permute(p,perms);
    Driver_for_Polyhedral_Continuation
      (file,permp,0,byebye,qq,qft,solsft,tosolve,ranstart,contrep);
   -- HANDLING THE UNMIXED AND THE MIXED CASE SEPARATELY :
    if not tosolve
     then if r = 1
           then Compute_Triangulation;
           else if caytrick
                 then Compute_Cayley_Triangulation;
                 else Compute_Mixed_Subdivision;
                end if;
          end if;
     else Solve_Coefficient_System;
    end if;
  end Driver_for_Dynamic_Mixed_Volume_Computation;

end Drivers_for_Dynamic_Lifting;
