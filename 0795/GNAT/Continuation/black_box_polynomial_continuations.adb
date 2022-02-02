with integer_io;                        use integer_io;
with File_Operations;                   use File_Operations;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Float_Vectors_of_Vectors;
with Complex_Numbers,Complex_Norms;     use Complex_Numbers,Complex_Norms;
with Complex_Vectors;
with Complex_Numbers_io;                use Complex_Numbers_io;
with Random_Number_Generators;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Solutions_io;                      use Solutions_io;
with Scaling;                           use Scaling;
with Projective_Transformations;        use Projective_Transformations;

with Continuation_Parameters;
with Continuation_Parameters_io;
with Homotopy,Process_io;               use Process_io;
with Increment_and_Fix_Continuation;    use Increment_and_Fix_Continuation;
with Root_Refiners;                     use Root_Refiners;
with Scanners_for_Continuation;         use Scanners_for_Continuation;

package body Black_Box_Polynomial_Continuations is

  procedure Scan_Input
                ( infile,outfile : in file_type; p,q : in out Link_to_Poly_Sys;
                  sols : in out Solution_List; arti : out boolean ) is

  -- DESCRIPTION :
  --   Scans the input file for target system and, if the homotopy is 
  --   artificial (in that case arti = true, otherwise arti = false),
  --   for a start system.  In both cases, start solutions are required.

    found,artificial : boolean;

  begin
    get(infile,p);
    put(outfile,p.all);
    artificial := (Number_of_Unknowns(p(p'first)) = p'last);
    if artificial
     then Scan_and_Skip(infile,"START SYSTEM",found);
          if found
           then get(infile,q);
                new_line(outfile);
                put_line(outfile,"THE START SYSTEM : ");
                new_line(outfile);
                put_line(outfile,q.all);
          end if;
    end if;
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found
     then get(infile,sols);
          new_line(outfile);
          put_line(outfile,"THE START SOLUTIONS : ");
          new_line(outfile);
          put(outfile,Length_Of(sols),Head_Of(sols).n,sols);
          new_line(outfile);
    end if;
    arti := artificial;
  end Scan_Input;

  procedure Set_Homotopy_Parameters
               ( file : in file_type; k : in out natural;
                 a,t : in out double_complex; prt : in out boolean ) is

  -- DESCRIPTION :
  --   Sets the default values for the homotopy parameters.

  begin
    k := 2;
    a := Random_Number_Generators.Random1;
    t := CMPLX(1.0);
    prt := false;
    new_line(file);
    put_line(file,"HOMOTOPY PARAMETERS :");
    put(file,"  k : "); put(file,k,2); new_line(file);
    put(file,"  a : "); put(file,a);   new_line(file);
    put(file,"  t : "); put(file,t);   new_line(file);
    if prt
     then put_line(file,"  projective transformation");
     else put_line(file,"  no projective transformation");
    end if;
  end Set_Homotopy_Parameters;

  procedure Tune_Continuation_Parameters ( infile,outfile : in file_type ) is

  -- DESCRIPTION :
  --   Scans the input file for continuation parameters and the
  --   output parameter.

  begin
    Continuation_Parameters.Tune(2);
    new_line(outfile);
    put_line(outfile,"****************** CURRENT CONTINUATION PARAMETERS "
      & "*****************");
    Continuation_Parameters_io.put(outfile);
    put_line(outfile,"***************************************************"
      & "*****************");
    Process_io.Set_Output_Code(nil);
  end Tune_Continuation_Parameters;

  procedure Black_Box_Refine
                  ( outfile : in file_type; p : in Poly_Sys;
                    artificial : in boolean; target : in double_complex;
                    k : in natural; sols,refsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Refines the roots in sols w.r.t. the system p.

    epsxa,epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    len,nb : natural := 0;

  begin
    if Length_Of(sols) > 0
     then 
       if artificial
        then if not Is_Null(sols) and then Head_Of(sols).n > p'last
              then Affine_Transformation(sols);
             end if;
             Reporting_Root_Refiner
               (outfile,p,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
        else declare
               pt : Poly_Sys(p'range);
             begin
               pt := Eval(p,target,k);
               Reporting_Root_Refiner
                 (outfile,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
               Clear(pt);
             end;
       end if;
    end if;
  end Black_Box_Refine;

  procedure Black_Box_Polynomial_Continuation
                 ( infile,outfile : in file_type; pocotime : out duration ) is

    p,q,sp : Link_to_Poly_Sys;
    sols,refsols : Solution_List;
    timer : timing_widget;
    k : natural := 0;
    a,target : double_complex;
    proj,artificial : boolean;
    rcond : double_float;
    scalecoeff : Complex_Vectors.Link_to_Vector;

    procedure Cont is
      new Reporting_Continue(Norm1,Homotopy.Eval,Homotopy.Diff,Homotopy.Diff);

   begin
    Scan_Input(infile,outfile,p,q,sols,artificial);
    scalecoeff := new Complex_Vectors.Vector(1..2*p'length);
    sp := new Poly_Sys(p'range);
    Copy(p.all,sp.all);
    Scale(sp.all,2,false,rcond,scalecoeff.all);
    Set_Homotopy_Parameters(outfile,k,a,target,proj);
    if artificial
     then Homotopy.Create(sp.all,q.all,k,a);
     else Homotopy.Create(sp.all,k);
          target := a;
    end if;
    Tune_Continuation_Parameters(infile,outfile);
    new_line(outfile);
    put_line(outfile,"THE SCALED SOLUTIONS :");
    new_line(outfile);
    tstart(timer);
    Cont(outfile,sols,proj,target);
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"continuation");
    pocotime := Elapsed_User_Time(timer);
    Scale(2,scalecoeff.all,sols);
    Clear(sp);
    Black_Box_Refine(outfile,p.all,artificial,target,k,sols,refsols);
  end Black_Box_Polynomial_Continuation;

  procedure Black_Box_Polynomial_Continuation
                  ( infile,outfile : in file_type;
                    p,q : in Poly_Sys; sols : in out Solution_List;
                    pocotime : out duration ) is

    refsols : Solution_List;
    timer : timing_widget;
    k : natural := 0;
    a,target : double_complex := CMPLX(0.0);
    proj : boolean;
    rcond : double_float;
    scalecoeff : Complex_Vectors.Vector(1..2*p'length);
    sp : Poly_Sys(p'range);

    procedure Cont is
      new Reporting_Continue(Norm1,Homotopy.Eval,Homotopy.Diff,Homotopy.Diff);

   begin
    Set_Homotopy_Parameters(outfile,k,a,target,proj);
    Copy(p,sp);
    Scale(sp,2,false,rcond,scalecoeff);
    Homotopy.Create(sp,q,k,a);
    Tune_Continuation_Parameters(infile,outfile);
    new_line(outfile);
    put_line(outfile,"THE SCALED SOLUTIONS :");
    new_line(outfile);
    tstart(timer);
    Cont(outfile,sols,proj,target);
    tstop(timer);
    new_line(outfile);
    print_times(outfile,timer,"continuation");
    pocotime := Elapsed_User_Time(timer);
    Scale(2,scalecoeff,sols);
    Clear(sp);
    Black_Box_Refine(outfile,p,true,target,k,sols,refsols);
    sols := refsols;
  end Black_Box_Polynomial_Continuation;

end Black_Box_Polynomial_Continuations;
