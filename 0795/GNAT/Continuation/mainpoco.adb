with Communications_with_User;            use Communications_with_User;
with text_io,integer_io,Numbers_io;       use text_io,integer_io,Numbers_io;
with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Complex_Numbers;                     use Complex_Numbers;
with Complex_Multivariate_Polynomials;    use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;          use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;       use Complex_Polynomial_Systems_io;
with Homotopy;
with Solutions,Solutions_io;              use Solutions,Solutions_io;
with Projective_Transformations;          use Projective_Transformations;
with Root_Refiners;                       use Root_Refiners;
with Drivers_for_Polynomial_Continuation;
 use Drivers_for_Polynomial_Continuation;
--with Bye_Bye_Message;

procedure mainpoco ( infilename,outfilename : in string ) is

  solsft,outft : file_type;
  lp : Link_to_Poly_Sys;
  sols,refsols : Solution_List;
  artificial,solsfile : boolean;
  k,len : natural;
  ans : character;
  tarre,tarim : double_float;
  target : double_complex;

  procedure Read_System ( filename : in string ) is
   
    file : file_type;

  begin
    if filename /= ""
     then Open(file,in_file,filename);
          get(file,lp);
          Close(file);
    end if;
  exception
    when others => 
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

begin
  Read_System(infilename);
  if lp = null
   then new_line; get(lp);
  end if;
  Create_Output_File(outft,outfilename);
  put(outft,lp.all); new_line(outft);
  new_line;
  put("Do you want the solutions on separate file ? (y/n) ");
  Ask_Yes_or_No(ans);
  if ans = 'y'
   then 
     put_line("Reading the name of the file to write the solutions on.");
     Read_Name_and_Create_File(solsft);
     solsfile := true;
   else
     solsfile := false;
  end if;
  artificial := (Number_of_Unknowns(lp(lp'first)) = lp'last);
  if artificial
   then Driver_for_Polynomial_Continuation(outft,lp.all,sols,target);
   else new_line;
        put("Give the index of the parameter : "); Read_Natural(k);
        new_line;
        put_line("Reading the target value of the continuation parameter.");
        put("Give the real part of the target : "); Read_Double_Float(tarre);
        put("Give the imaginary part of the target : ");
        Read_Double_Float(tarim);
        target := CMPLX(tarre,tarim);
        Driver_for_Polynomial_Continuation(outft,lp.all,k,target,sols);
  end if;
  if Length_Of(sols) > 0
   then declare
          epsxa,epsfa,tolsing : constant double_float := 10.0**(-8);
          nb : natural := 0;
        begin
          if artificial
           then
             if not Is_Null(sols) and then Head_Of(sols).n > lp'last
              then Affine_Transformation(sols);
             end if;
             if target = CMPLX(1.0)
              then
                if solsfile
                 then Reporting_Root_Refiner
                        (outft,lp.all,sols,refsols,epsxa,epsfa,tolsing,
                         nb,5,false);
                 else Reporting_Root_Refiner
                        (outft,lp.all,sols,epsxa,epsfa,tolsing,nb,5,false);
                end if;
              else
                declare
                  pt : Poly_Sys(lp'range);
                begin
                  pt := Homotopy.Eval(target);
                  if solsfile
                   then Reporting_Root_Refiner
                          (outft,pt,sols,refsols,epsxa,epsfa,tolsing,
                           nb,5,false);
                   else Reporting_Root_Refiner
                          (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,false);
                  end if;
                  Clear(pt);
                end;
             end if;
           else
             declare
               pt : Poly_Sys(lp'range);
             begin
               pt := Eval(lp.all,target,k);
               if solsfile
                then Reporting_Root_Refiner
                       (outft,pt,sols,refsols,epsxa,epsfa,tolsing,nb,5,false);
                else Reporting_Root_Refiner
                       (outft,pt,sols,epsxa,epsfa,tolsing,nb,5,false);
               end if;
               Clear(pt);
             end;
          end if;
        end;
  end if;
 -- put(outft,Bye_Bye_Message);
  Close(outft);
  if solsfile
   then len := Length_Of(refsols);
        if len > 0
         then put(solsft,len,Head_Of(refsols).n,refsols);
        end if;
        Close(solsft);
  end if;
end mainpoco;
