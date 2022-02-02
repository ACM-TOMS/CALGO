with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;

with text_io,integer_io;                use text_io,integer_io;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers;                   use Complex_Numbers;
with Integer_Vectors;                   use Integer_Vectors;
with Integer_Vectors_io;                use Integer_Vectors_io;
with Float_Vectors,Float_Vectors_io;    use Float_Vectors_io;
with Complex_Vectors;
with Complex_Vectors_io;                use Complex_Vectors_io;
with Integer_Vectors_of_Vectors;
with Float_Vectors_of_Vectors;
with Complex_Vectors_of_Vectors;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;       use Lists_of_Integer_Vectors_io;

with Integer_Vectors_of_Vectors_io;     use Integer_Vectors_of_Vectors_io;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Float_Vector_Lists;      use Arrays_of_Float_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;

with Power_Lists;                       use Power_Lists;
with Integer_Faces_of_Polytope;
with Float_Faces_of_Polytope;
with Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;     use Integer_Mixed_Subdivisions_io;
with Float_Mixed_Subdivisions;
with Float_Mixed_Subdivisions_io;       use Float_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;          use Mixed_Volume_Computation;

with Float_Integer_Convertors;          use Float_Integer_Convertors;
with Integer_Lifting_Functions;         use Integer_Lifting_Functions;
with Integer_Lifting_Utilities;         use Integer_Lifting_Utilities;
with Float_Lifting_Functions;           use Float_Lifting_Functions;
with Float_Lifting_Utilities;           use Float_Lifting_Utilities;
with Integer_Pruning_Methods;           use Integer_Pruning_Methods;
with Float_Pruning_Methods;             use Float_Pruning_Methods;

with Solutions,Solutions_io;            use Solutions,Solutions_io;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;
with Complex_Multivariate_Laurent_Polynomials;
 use Complex_Multivariate_Laurent_Polynomials;
with Exponent_Vectors;                  use Exponent_Vectors;
with Polynomial_to_Laurent_Converters;  use Polynomial_to_Laurent_Converters;
with Laurent_to_Polynomial_Converters;  use Laurent_to_Polynomial_Converters;
with Laurent_Jacobi_Matrices;           use Laurent_Jacobi_Matrices;

with Driver_for_Criterion;
with Drivers_for_Lifting_Functions;     use Drivers_for_Lifting_Functions;
with Pruning_Statistics;
with Integer_Polyhedral_Continuation;   use Integer_Polyhedral_Continuation;
with Float_Polyhedral_Continuation;     use Float_Polyhedral_Continuation;
with Driver_for_Polyhedral_Continuation;

package body Drivers_for_Static_Lifting is

  procedure Static_Lifting_Info is

    i : array(1..11) of string(1..65);

  begin
    i( 1):="  Static lifting is a  general  procedure  to  construct  regular";
    i( 2):="mixed  subdivisions  of  a tuple of polytopes.   For mixed-volume";
    i( 3):="computation, only those cells that are  spanned  by  a  tuple  of";
    i( 4):="edges  contribute  to  the mixed volume.  These cells are the so-";
    i( 5):="called mixed cells in the subdivision.  The collection  of  mixed";
    i( 6):="cells  is  computed  efficiently by pruning in the tree of lifted";
    i( 7):="edge-edge combinations.                                          ";
    i( 8):="  These mixed cells provide the start systems in  the  polyhedral";
    i( 9):="homotopy methods used to solve a random coefficient start system.";
    i(10):="Recursion is applied in case the lifting does not induce at  once";
    i(11):="a fine mixed subdivision.                                        ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Static_Lifting_Info;

  procedure Driver_for_Mixed_Volume_Computation 
               ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                 q : out Poly_Sys; qsols : out Solution_List;
                 mv : out natural ) is

    welcome : constant string := "Mixed-Volume Computation by Static Lifting";

    permp : Poly_Sys(p'range);          -- permuted system
    tmv,nbcells : natural;
    mix : Integer_Vectors.Link_to_Vector;
    imixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    fmixsub : Float_Mixed_Subdivisions.Mixed_Subdivision;

   -- the switches, to do or not to do, if true then :

    compmix   : boolean;   -- compute the type of mixture
    compmisu  : boolean;   -- compute the subdivision
    report    : boolean;   -- report during creation of the subdivision
    misufile  : boolean;   -- put the subdivision on separate file 
    tosolve   : boolean;   -- solve the system 
    ranstart  : boolean;   -- construct random coefficient start system
    contrep   : boolean;   -- intermediate output during continuation

   -- the files :

    solsft,gft : file_type;

  -- AUXILIAIRY :

    function Select_Terms
                 ( p : Poly_Sys; mix : Integer_Vectors.Vector;
                   pts : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                 return Poly_Sys is

    -- DESCRIPTION :
    --   Returns those terms of p whose exponents belong to the lists pts.
    --   The type of mixture is in mix and pts'range = mix'range.

      res : Poly_Sys(p'range);
      ind : natural := 0;

    begin
      for i in pts'range loop
        for j in 1..mix(i) loop
          ind := ind + 1;
          res(ind) := Select_Terms(p(ind),pts(i));
        end loop;
      end loop;
      return res;
    end Select_Terms;

  -- REPORT DURING CREATION :

    procedure Write_Cell ( mic : in Integer_Mixed_Subdivisions.Mixed_Cell;
                           continue : out boolean ) is

      vol : natural;

    begin
      nbcells := nbcells + 1;
      put(file,"Cell no. "); put(file,nbcells,1); put_line(file," : ");
      put(file," normal to cell : "); put(file,mic.nor); new_line(file);
      put_line(file," the points in the cell : ");
      for k in mic.pts'range loop
        put(file,"  component "); put(file,k,1); put(file," with ");
        put(file,Length_Of(mic.pts(k)),1); put_line(file," points :");
        put(file,mic.pts(k));
      end loop;
      vol := Mixed_Volume(p'length,mix.all,mic);
      put(file," with volume addition : ");
      put(file,tmv,1); put(file," + "); put(file,vol,1);
      tmv := tmv + vol; put(file," = "); put(file,tmv,1); new_line(file);
      continue := true;
    end Write_Cell;
    procedure Report_and_Create1_CS is new Gen1_Create_CS(Write_Cell);

    procedure Compute_Mixture 
              ( file : in file_type; n : in natural;
                compmixture : in boolean;
                supports : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in out Integer_Vectors.Link_to_Vector;
                permp : out Poly_Sys ) is

    -- DESCRIPTION :
    --   Computes the type of mixture and permutes if necessary,
    --   the equations in the polynomial system.

      perm : Integer_Vectors.Link_to_Vector;

    begin
      if compmixture
       then Compute_Mixture(supports,mix,perm);
       else perm := Compute_Permutation(n,mix.all,supports);
      end if;
      permp := Permute(p,perm);  Clear(perm);
      new_line(file);
      put(file,"TYPE OF MIXTURE : "); put(file,"#supports : ");
      put(file,mix'last,1);
      put(file," occurrences : "); put(file,mix);
      new_line(file);
    end Compute_Mixture;

    function Expand
                ( mix : Vector;
                  points : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    -- DESCRIPTION :
    --   Returns a tuple of expanded lists, according to the type of mixture.

      sum : integer := Integer_Vectors.Sum(mix);
      res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..sum);
      cnt : natural := 0;

    begin
      for i in mix'range loop
        for j in 1..mix(i) loop
          cnt := cnt + 1;
          res(cnt) := points(i);
        end loop;
      end loop;
      return res;
    end Expand;

    procedure Write_Cardinalities
               ( file : in file_type; mix,card : in Integer_Vectors.Vector ) is
 
    begin
      new_line(file);
      put_line(file,"CARDINALITIES OF THE LIFTED FACES :");
      new_line(file);
      for i in card'range loop
        put(file,"  #lifted "); put(file,mix(i),1);
        put(file,"-faces of polytope "); put(file,i,1);
        put(file," : "); put(file,card(i),1); new_line(file);
      end loop;
    end Write_Cardinalities;

    procedure Create_Mixed_Cells
             ( file : in file_type; n : in natural; mix : in Vector;
               mixpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lifted : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   The pruning algorithm will be applied to compute the mixed cells.

      use Integer_Faces_of_Polytope;
      afa : Array_of_Faces(mix'range);
      cardafa : Integer_Vectors.Vector(mix'range);
      nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);
      timer : timing_widget;

    begin
      tstart(timer);
      for i in afa'range loop
        afa(i) := Create_Lower(mix(i),n+1,lifted(i));
      end loop;
      tstop(timer);
      for i in afa'range loop
        cardafa(i) := Length_Of(afa(i));
      end loop;
      Write_Cardinalities(file,mix,cardafa);
      new_line(file);
      print_times(file,timer,"Creation of the faces of lower hull");
      new_line(file);
      put_line(file,"PRUNING TO EXTRACT THE MIXED CELLS :");
      tstart(timer);
      if report
       then nbcells := 0; tmv := 0;
            Report_and_Create1_CS(n,mix,afa,lifted,nbsucc,nbfail,mixsub);
       else Create_CS(n,mix,afa,lifted,nbsucc,nbfail,mixsub);
      end if;
      tstop(timer);
      Pruning_Statistics(file,nbsucc,nbfail);
      new_line(file);
      print_times(file,timer,"Pruning for Mixed Cells");
    end Create_Mixed_Cells;

    procedure Create_Mixed_Cells
              ( file : in file_type;
                n : in natural; mix : in Integer_Vectors.Vector;
                fltsup : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
                lilifu : in Float_Vectors_of_Vectors.Link_to_Vector;
                lifsup : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
                fltsub : in out Float_Mixed_Subdivisions.Mixed_Subdivision ) is

      use Float_Vectors_of_Vectors;
      use Float_Faces_of_Polytope;
      timer : timing_widget;
      tol : constant double_float := 10.0**(-8);
      supfa,liffaces : Array_of_Faces(mix'range);
      cardafa : Integer_Vectors.Vector(mix'range);
      nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);

    begin
      tstart(timer);
      if lilifu /= null
       then for i in supfa'range loop
              supfa(i) := Create(mix(i),n,fltsup(i),tol);
              liffaces(i) := Linear_Lift(supfa(i),lilifu(i).all);
            end loop;
       else for i in liffaces'range loop
              liffaces(i) := Create_Lower(mix(i),n+1,lifsup(i),tol);
            end loop;
      end if;
      tstop(timer);
      for i in liffaces'range loop
        cardafa(i) := Length_Of(liffaces(i));
      end loop;
      Write_Cardinalities(file,mix,cardafa);
      new_line(file);
      print_times(file,timer,"Creation of the faces of lower hull");
      new_line(file);
      put_line(file,"PRUNING TO EXTRACT THE MIXED CELLS :");
      tstart(timer);
      Create(n,mix,liffaces,lifsup,tol,nbsucc,nbfail,fltsub);
      tstop(timer);
      Pruning_Statistics(file,nbsucc,nbfail);
      new_line(file);
      print_times(file,timer,"Pruning for Mixed Cells.");
    end Create_Mixed_Cells;

  procedure Polyhedral_Homotopy_Continuation
                 ( file : in file_type; qq : in out Poly_Sys;
                 qqsols : in out Solution_List;
                 lifted : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 mixsub : in Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   This procedure constructs and solves a random coefficient start system.
    --   Only those monomials whose exponents occur in the subdivision are 
    --   contained in the start system.

      timer : timing_widget;
      n : constant natural := p'length;
      lifted_lq,lq : Laur_Sys(q'range);
      h : Eval_Coeff_Laur_Sys(q'range);
      c : Complex_Vectors_of_Vectors.Vector(h'range);
      e : Exponent_Vectors_Array(h'range);
      j : Eval_Coeff_Jacobi(h'range,h'first..h'last+1);
      m : Mult_Factors(j'range(1),j'range(2));

    begin
      new_line(file);
      put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION :");
      lq := Polynomial_to_Laurent_System(qq);
      lifted_lq := Perform_Lifting(n,mix.all,lifted,lq);
      Clear(lq);
      lq := Eval(lifted_lq,CMPLX(1.0),n+1);
      qq := Laurent_to_Polynomial_System(lq);
      h := Create(lq);
      for i in c'range loop
        declare
          coeff_lq : constant Complex_Vectors.Vector := Coeff(lq(i));
        begin
          c(i) := new Complex_Vectors.Vector(coeff_lq'range);
          for k in coeff_lq'range loop
            c(i)(k) := coeff_lq(k);
          end loop;
        end;
      end loop;
      e := Create(lq);
      Create(lq,j,m);
      tstart(timer);
      if contrep
       then --Mixed_Solve(file,lifted_lq,mix.all,mixsub,qqsols);
            Mixed_Solve(file,lifted_lq,lifted,h,c,e,j,m,mix.all,mixsub,qqsols);
       else --Mixed_Solve(lifted_lq,mix.all,mixsub,qqsols);
            Mixed_Solve(lifted_lq,lifted,h,c,e,j,m,mix.all,mixsub,qqsols);
      end if;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Polyhedral Continuation");
      q := qq; qsols := qqsols;
  end Polyhedral_Homotopy_Continuation;

    procedure Polyhedral_Homotopy_Continuation
                ( file : in file_type;
                  n : in natural; mix : in Integer_Vectors.Vector;
                  qq : in Poly_Sys; qqsols : in out Solution_List;
                  lifsup : in Arrays_of_Float_Vector_Lists.Array_of_Lists;
                  fltsub : in Float_Mixed_Subdivisions.Mixed_Subdivision;
                  contrep : in boolean ) is

      lq : Laur_Sys(qq'range) := Polynomial_to_Laurent_System(qq);
      h : Eval_Coeff_Laur_Sys(qq'range);
      c : Complex_Vectors_of_Vectors.Vector(h'range);
      e : Exponent_Vectors_Array(h'range);
      j : Eval_Coeff_Jacobi(h'range,h'first..h'last+1);
      m : Mult_Factors(j'range(1),j'range(2));
      timer : timing_widget;

    begin
      new_line(file);
      put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION :");
      h := Create(lq);
      for i in c'range loop
        declare
          coeff_lq : constant Complex_Vectors.Vector := Coeff(lq(i));
        begin
          c(i) := new Complex_Vectors.Vector(coeff_lq'range);
          for k in coeff_lq'range loop
            c(i)(k) := coeff_lq(k);
          end loop;
        end;
      end loop;
      e := Create(lq);
      Create(lq,j,m);
      tstart(timer);
      if contrep
       then Mixed_Solve(file,lq,lifsup,h,c,e,j,m,mix,fltsub,qqsols);
       else Mixed_Solve(lq,lifsup,h,c,e,j,m,mix,fltsub,qqsols);
      end if;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"Polyhedral Homotopy Continuation");
      q := qq; qsols := qqsols;
    end Polyhedral_Homotopy_Continuation;

    procedure Volume_Computation
             ( file : in file_type;
               n : in natural; mix : in Vector;
               lifpts : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

    -- DESCRIPTION :
    --   Computes the volumes of the mixed cells in mixsub.

      timer : timing_widget;
      mixvol : natural;

    begin
      new_line(file);
      put_line(file,"VOLUMES OF MIXED CELLS :");
      new_line(file);
      tstart(timer);
      put(file,n,mix,mixsub,mixvol);
      tstop(timer);
      put(file,"The total mixed volume equals ");
      put(file,mixvol,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Volume computation of mixed cells");
      if compmisu
       then mixsub := Integer_Mixed_Subdivisions.Create(lifpts,mixsub);
            new_line(file);
            put_line(file,"VOLUMES OF MIXED CELLS, AFTER REFINEMENT :");
            new_line(file);
            tstart(timer);
            put(file,n,mix,mixsub,mixvol);
            tstop(timer);
            put(file,"The total mixed volume equals ");
            put(file,mixvol,1); new_line(file);
            new_line(file);
            print_times(file,timer,"Volume computation of mixed cells");
      end if;
      mv := mixvol;
    end Volume_Computation;

    procedure Volume_Computation
               ( file : in file_type; n : in natural; mix : in Vector;
                 mixsub : in out Float_Mixed_Subdivisions.Mixed_Subdivision ) is

      timer : timing_widget;
      mixvol : natural;

    begin
      new_line(file);
      put_line(file,"THE MIXED SUBDIVISION : ");
      new_line(file);
      tstart(timer);
      Float_Mixed_Subdivisions_io.put(file,n,mix,mixsub,mixvol);
      tstop(timer);
      put(file,"The mixed volume equals : ");
        put(file,mixvol,1); new_line(file);
      new_line(file);
      print_times(file,timer,"Volume computation of mixed cells");
      mv := mixvol;
    end Volume_Computation;

    procedure Data_Management
                    ( mix : in out Integer_Vectors.Link_to_Vector;
                      compmisu,compmix,fltlif : out boolean ) is

    -- DESCRIPTION :
    --   This procedure allows to use previously computed mixed subdivisions.

      ans : character;
      m : natural;

    begin
      new_line;
      put("Do you already have a mixed subdivision ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put("Induced by integer or floating-point lifting (i/f) ");
            Ask_Alternative(ans,"if");
            fltlif := (ans = 'f');
            declare
              insubft : file_type;
              nn : natural;
            begin
              put_line("Reading the name of the input file.");
              Read_Name_and_Open_File(insubft);
              if ans = 'f'
               then get(insubft,nn,m,mix,fmixsub);
               else get(insubft,nn,m,mix,imixsub);
              end if;
              Close(insubft);
              new_line(file);
              put_line(file,"Mixed subdivision supplied by user.");
              new_line(file);
              compmisu := false; compmix := false;
            exception
              when DATA_ERROR
                => put_line("Data not in correct format.  Will ignore it...");
                   Close(insubft);
            end;
       else compmisu := true;
            put("Do you want to enforce a type mixture ? (y/n) ");
            Ask_Yes_or_No(ans);
            if ans = 'y'
             then put("Give number of different supports : "); get(m);
                  put("Give vector indicating occurrences : ");
                  get(m,mix);
                  compmix := false;
             else compmix := true;
            end if;
      end if;
    end Data_Management;

    procedure Main_Driver is

      n : constant natural := p'length;
      sp,qq : Poly_Sys(p'range);
      qqsols : Solution_List;
      totaltimer : timing_widget;
      points : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
      ilili : Integer_Vectors_of_Vectors.Link_to_Vector;
      flili : Float_Vectors_of_Vectors.Link_to_Vector;
      fltlif : boolean;
      mixpts,mixpts1,ilifpts 
             : Arrays_of_Integer_Vector_Lists.Link_to_Array_of_Lists;
      fpts,flifpts : Arrays_of_Float_Vector_Lists.Link_to_Array_of_Lists;
      ans : character;

      outsubft : file_type;

      procedure Driver_for_Integer_Lifting is
      begin
        if compmisu
         then new_line;
              put("Do you want intermediate output on file,");
              put(" during creation ? (y/n) ");
              Ask_Yes_or_No(ans);
              report := (ans = 'y');
              put("Do you want the mixed cells on separate file ? (y/n) ");
              Ask_Yes_or_No(ans);
              misufile := (ans = 'y');
              if misufile
               then put_line("Reading the name of the output file.");
                    Read_Name_and_Create_File(outsubft);
              end if;
         else ilifpts.all := Induced_Lifting(n,mix.all,points,imixsub);
        end if;
        Driver_for_Polyhedral_Continuation
          (file,sp,0,byebye,qq,gft,solsft,tosolve,ranstart,contrep);
        if compmisu
         then Create_Mixed_Cells(file,n,mix.all,mixpts.all,ilifpts.all,imixsub);
        end if;
        if not Integer_Mixed_Subdivisions.Is_Null(imixsub)
         then Volume_Computation(file,n,mix.all,ilifpts.all,imixsub);
              if compmisu and then misufile
               then put(outsubft,n,mix.all,imixsub);
              end if;
              if tosolve
               then Polyhedral_Homotopy_Continuation
                      (file,qq,qqsols,ilifpts.all,imixsub);
              end if;
        end if;
      end Driver_for_Integer_Lifting;

      procedure Driver_for_Float_Lifting is
      begin
        if compmisu
         then new_line;
              put("Do you want the mixed cells on separate file ? (y/n) ");
              Ask_Yes_or_No(ans);
              misufile := (ans = 'y');
              if misufile
               then put_line("Reading the name of the output file.");
                    Read_Name_and_Create_File(outsubft);
              end if;
         else flifpts.all := Induced_Lifting(n,mix.all,fpts.all,fmixsub);
        end if;
        Driver_for_Polyhedral_Continuation
          (file,sp,0,byebye,qq,gft,solsft,tosolve,ranstart,contrep);
        if compmisu
         then Create_Mixed_Cells
                 (file,n,mix.all,fpts.all,flili,flifpts.all,fmixsub);
              if misufile
               then put(outsubft,n,mix.all,fmixsub);
              end if;
        end if;
        if not Float_Mixed_Subdivisions.Is_Null(fmixsub)
         then Volume_Computation(file,n,mix.all,fmixsub);
              if tosolve
               then Polyhedral_Homotopy_Continuation
                      (file,n,mix.all,qq,qqsols,flifpts.all,fmixsub,contrep);
              end if;
        end if;
      end Driver_for_Float_Lifting;

    begin
      new_line; put_line(welcome);
      tstart(totaltimer);
  
      points := Construct_Power_Lists(p);

      Data_Management(mix,compmisu,compmix,fltlif);

      if compmisu
       then
         Compute_Mixture(file,n,compmix,points,mix,permp);
         mixpts1 := new Arrays_of_Integer_Vector_Lists.
                        Array_of_Lists'(Typed_Lists(mix.all,points));
  
         Driver_for_Criterion(file,mixpts1.all);

         mixpts := new Arrays_of_Integer_Vector_Lists.
                       Array_of_Lists'(Expand(mix.all,mixpts1.all));
         Compute_Mixture(file,n,compmix,mixpts.all,mix,permp);
         mixpts := new Arrays_of_Integer_Vector_Lists.
                       Array_of_Lists'(Typed_Lists(mix.all,mixpts.all));

         ilifpts := new Arrays_of_Integer_Vector_Lists.
                        Array_of_Lists(mix'range);
         fpts := new Arrays_of_Float_Vector_Lists.Array_of_Lists(mix'range);
         flifpts := new Arrays_of_Float_Vector_Lists.Array_of_Lists(mix'range);

         sp := Select_Terms(permp,mix.all,mixpts.all);

         new_line;
         Driver_for_Lifting_Functions
           (file,sp,mixpts.all,fltlif,fpts.all,ilifpts.all,
                 flifpts.all,ilili,flili);
       else
         mixpts := new Arrays_of_Integer_Vector_Lists.
                       Array_of_Lists'(Typed_Lists(mix.all,points));
         ilifpts := new Arrays_of_Integer_Vector_Lists.
                        Array_of_Lists(mix'range);
         fpts := new Arrays_of_Float_Vector_Lists.Array_of_Lists(points'range);
         fpts.all := Convert(points);
         flifpts := new Arrays_of_Float_Vector_Lists.Array_of_Lists(mix'range);
         sp := p;
      end if;
      if fltlif
       then Driver_for_Float_Lifting;
       else Driver_for_Integer_Lifting;
      end if;

      if tosolve
       then new_line(file);
            put_line(file,"THE RANDOM COEFFICIENT START SYSTEM :");
            new_line(file);
            put_line(file,qq);
            new_line(file);
            put_line(file,"THE START SOLUTIONS :");
            new_line(file);
            put(file,Length_Of(qqsols),qq'length,qqsols);
            if ranstart
             then new_line(gft);
                  put_line(gft,"THE SOLUTIONS : "); new_line(gft);
                  put(gft,Length_Of(qqsols),n,qqsols);
                  Close(gft);
             else put(solsft,Length_Of(qqsols),n,qqsols);
                  Close(solsft);
            end if;
       end if;
 
       tstop(totaltimer);
       new_line(file);
       print_times(file,totaltimer,"All Computations");
     end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Mixed_Volume_Computation;

end Drivers_for_Static_Lifting;
