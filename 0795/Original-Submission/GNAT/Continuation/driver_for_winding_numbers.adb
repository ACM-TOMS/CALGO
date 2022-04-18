with integer_io,Numbers_io;              use integer_io,Numbers_io;
with Communications_with_User;           use Communications_with_User;
with File_Operations;                    use File_Operations;
with Timing_Package;                     use Timing_Package;

with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Complex_Vectors,Complex_Norms;      use Complex_Vectors,Complex_Norms;
with Homotopy;
with Solutions_io;                       use Solutions_io;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems_io;      use Complex_Polynomial_Systems_io;

with Continuation_Parameters;            use Continuation_Parameters;
with Continuation_Data;                  use Continuation_Data;
with Path_Trackers;                      use Path_Trackers;

with Driver_for_Homotopy_Construction;
with Drivers_for_Polynomial_Continuation;
 use Drivers_for_Polynomial_Continuation;

procedure Driver_for_Winding_Numbers
             ( file : in file_type; p : in Poly_Sys;
               sols : in out Solution_List ) is

  use Floating_Point_Numbers.double_float_io;

  infile : file_type;
  timer : timing_widget;
  pp,q : Poly_Sys(p'range);
  found,proj : boolean := false;
  target : double_complex;
  ans : character;
  oc,max_wc : natural;

  procedure Write_Statistics_and_Condition
                ( file : in file_type; i,nstep,nfail,niter,nsyst : in natural;
                  rcond : in double_float ) is

  -- DESCRIPTION :
  --   Writes the computing statistics of the ith path on file, followed
  --   by the estimate for the inverse of the condition number.

  begin
    put(file,"========================================");
    put_line(file,"===================================");
    put(file,"== "); put(file,i,1); put(file," = ");
    put(file," #step : "); put(file,nstep,3);
    put(file," #fail : "); put(file,nfail,2);
    put(file," #iter : "); put(file,niter,3);
    if nsyst /= niter
     then put(file," #syst : "); put(file,nsyst,3);
    end if;
    put(file," = ");
    put(file," rco : "); put(file,rcond,2,3,3);
    put_line(file," ==");
  end Write_Statistics_and_Condition;

  procedure Winding_Numbers
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean; target : in double_complex;
                 mwc : in natural ) is

  -- DESCRIPTION :
  --   Computes the winding number for the given list of solutions.

  -- REQUIRED :
  --   The homotopy is already defined and stored in the package homotopy.

    sa : Solu_Info_Array(1..Length_Of(sols)) := Deep_Create(sols);
    pen : Pred_Pars := Continuation_Parameters.Create_End_Game;
    cen : Corr_Pars := Continuation_Parameters.Create_End_Game;
    tol : constant double_float := 10.0**(-10);
    epslop : constant double_float := 10.0**(-6);
    wc : natural;
    sum,allsum : Complex_Vectors.Vector(sa(sa'first).sol.v'range);

    procedure CCont2 is
      new Circular_Single_Conditioned_Reporting_Continue
            (Norm1,Homotopy.Eval,Homotopy.Diff,Homotopy.Diff);

  begin
    for i in sa'range loop
      if modulus(CMPLX(1.0) - sa(i).sol.t) > epslop
       then
         CCont2(file,sa(i),target,tol,epslop,wc,mwc,sum,allsum,false,pen,cen);
         sa(i).sol.m := wc;
         sa(i).sol.v := allsum;
         sa(i).sol.t := target;
         Write_Statistics_and_Condition
           (file,i,sa(i).nstep,sa(i).nfail,sa(i).niter,sa(i).nsyst,sa(i).rcond);
         put(file,sa(i).sol.all);
      end if;
    end loop;
    Clear(sols);
    sols := Shallow_Create(sa);
  end Winding_Numbers;

begin
 -- READING MAXIMUM WINDING NUMBER :
  new_line;
  put("Give the maximum winding number : "); Read_Natural(max_wc);
 -- READING THE START SYSTEM :
  new_line;
  put_line("Reading the name of the file for start system.");
  Read_Name_and_Open_File(infile); get(infile,q);
  new_line;
  put_line(file,"THE START SYSTEM :");
  new_line(file); put(file,q); new_line(file);
 -- CONSTRUCTING THE HOMOTOPY AND TUNING OF PARAMETERS :
  Copy(p,pp);
  Driver_for_Homotopy_Construction(file,pp,q,sols,target);
  proj := (Number_of_Unknowns(q(q'first)) > q'last);
  Driver_for_Continuation_Parameters(file);
  new_line;
  put("Do you want intermediate output during continuation ? (y/n) ");
  Ask_Yes_or_No(ans);
  if ans = 'y'
   then Driver_for_Process_io(file,oc);
  end if;
 -- COMPUTATION OF WINDING NUMBERS :
  new_line;
  put_line("No more input desired.  Computing winding numbers ...");
  put_line("The program can now run in the background.");
  new_line;
  tstart(timer);
  Winding_Numbers(file,sols,proj,target,max_wc);
  tstop(timer);
  new_line(file);
  print_times(file,timer,"computation of winding numbers");
  new_line(file);
end Driver_for_Winding_Numbers;
