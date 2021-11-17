with integer_io;                          use integer_io;
with Timing_Package;                      use Timing_Package;
with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Root_Refiners,Scaling;               use Root_Refiners,Scaling;
with Projective_Transformations;          use Projective_Transformations;

procedure Driver_for_Root_Refining
             ( file : in file_type; scalp,p : in Poly_Sys; basis : in natural;
               scalvec : in Link_to_Vector; sols : in out Solution_List ) is

  numb : natural;
  epsxa,epsfa : constant double_float := 10.0**(-8);
  tolsing : constant double_float := 10.0**(-8);
  timer : timing_widget;
  len : constant natural := Length_Of(sols);

begin
  if (len /= 0) and then Head_Of(sols).n > p'last
   then Affine_Transformation(sols);
  end if;
  if scalvec /= null
   then put_line(file,"ROOT REFINING ON THE SCALED SYSTEM :");
        tstart(timer);
        numb := 0;
        Reporting_Root_Refiner
          (file,scalp,sols,epsxa,epsfa,tolsing,numb,5,false);
        tstop(timer);
        new_line(file);
        print_times(file,timer,"Root Refining on the Scaled System");
        Scale(basis,scalvec.all,sols);
  end if;
  tstart(timer);
  numb := 0;
  Reporting_Root_Refiner(file,p,sols,epsxa,epsfa,tolsing,numb,5,false);
  tstop(timer);
  new_line(file);
  print_times(file,timer,"Root Refining");
end Driver_for_Root_Refining;
