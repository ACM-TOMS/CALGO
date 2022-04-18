with integer_io,Numbers_io;               use integer_io,Numbers_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Float_Vectors,Complex_Vectors;
with Complex_Norms,Complex_Matrices;      use Complex_Norms,Complex_Matrices;
with Substitutors;                        use Substitutors;
with Complex_Polynomial_Systems;          use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;       use Complex_Polynomial_Systems_io;
with Homotopy;
with Projective_Transformations;          use Projective_Transformations;
--with Toric_Compactification_Data;
with Continuation_Data;                   use Continuation_Data;
with Continuation_Parameters;             use Continuation_Parameters;
with Increment_and_Fix_Continuation;      use Increment_and_Fix_Continuation;

package body Drivers_for_Path_Directions is

-- AUXILIARIES TO INSTANTIATE :

  function Scale ( v : Float_Vectors.Vector ) return Float_Vectors.Vector is

  -- DESCRIPTION :
  --   Divides the vector by its largest component.

    res : Float_Vectors.Vector(v'range);
    tol : constant double_float := 10.0**(-12);
    ind : integer := v'first;
    max : double_float := abs(v(ind));

  begin
    for i in v'range loop
      if abs(v(i)) > max
       then max := abs(v(i));
            ind := i;
      end if;
    end loop;
    if max > tol
     then for i in v'range loop
            res(i) := v(i)/max;
          end loop;
    end if;
    return res;
  end Scale;

  function Toric_Evaluator ( x : Complex_Vectors.Vector; t : double_complex )
                           return Complex_Vectors.Vector is
  begin
   -- if not Toric_Compactification_Data.Numeric_Empty
   --  then return Toric_Compactification_Data.Eval(x,t);
   --  else
         return Homotopy.Eval(x,t);
   -- end if;
  end Toric_Evaluator;

  function Toric_Differentiator
               ( x : Complex_Vectors.Vector; t : double_complex )
               return Complex_Vectors.Vector is
  begin
    --if not Toric_Compactification_Data.Numeric_Empty
    -- then return Toric_Compactification_Data.Eval(x,t);
    -- else
         return Homotopy.Eval(x,t);
    --end if;
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : Complex_Vectors.Vector; t : double_complex )
               return Complex_Matrices.Matrix is
  begin
    --if not Toric_Compactification_Data.Numeric_Empty
    -- then return Toric_Compactification_Data.Diff(x,t);
    -- else
         return Homotopy.Diff(x,t);
    --end if;
  end Toric_Differentiator;

-- TARGET ROUTINES :

  procedure Init_Path_Directions
               ( n,nv : in natural;
                 v : in out Float_Vectors_of_Vectors.Link_to_Vector;
                 errv : in out Float_Vectors.Link_to_Vector ) is

  begin
    v := new Float_Vectors_of_Vectors.Vector(1..nv);
    for i in v'range loop
      v(i) := new Float_Vectors.Vector'(1..n => 0.0);
    end loop;
    errv := new Float_Vectors.Vector'(1..nv => 1.0);
  end Init_Path_Directions;

  procedure Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj,report : in boolean;
                 v : in out Float_Vectors_of_Vectors.Vector;
                 errv : in out Float_Vectors.Vector;
                 target : in double_complex ) is

  -- DESCRIPTION :
  --   Performs the continuation with online toric compactifications.

    timer : timing_widget;
    h : constant Poly_Sys := Homotopy.Eval(target); --Homotopy.Homotopy_System;
    n : constant natural := h'length;
    hh : Poly_Sys(h'range);

    procedure Sil_Cont is
      new Silent_Toric_Continue(Norm1,Toric_Evaluator,
                                Toric_Differentiator,Toric_Differentiator);
    procedure Rep_Cont is
      new Reporting_Toric_Continue(Norm1,Toric_Evaluator,
                                   Toric_Differentiator,Toric_Differentiator);
  begin
  --  if proj
  --   then Copy(h,hh); Toric_Compactification_Data.Initialize(hh,n);
  --   else Toric_Compactification_Data.Symbolic_Initialize(h,n+1);
  --  end if;
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,v,errv,target);
     else Sil_Cont(sols,proj,v,errv,target);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"toric continuation"); new_line(file);
  end Toric_Continue;

  procedure Write_Directions 
               ( file : in file_type;
                 v : in Float_Vectors_of_Vectors.Vector;
                 errv : in Float_Vectors.Vector ) is

    use Floating_Point_Numbers.double_float_io;
    use Float_Vectors;

    procedure Write ( file : in file_type; v : in Vector ) is
    begin
      for i in v'range loop
        put(file,v(i)); new_line(file);
      end loop;
    end Write;

    procedure Write_Direction
                 ( file : in file_type;
                   v : in Vector; error : in double_float; i : natural ) is

     -- sv : Vector(v'range) := Scale(v);

    begin
      put(file,"Computed direction of path ");
      put(file,i,1); put_line(file," :"); Write(file,v);
      put(file,"with error : "); put(file,error); new_line(file);
     -- put(file,"Scaled direction of path ");
     -- put(file,i,1); put_line(file," :"); Write(file,sv);
     -- if not Toric_Compactification_Data.Symbolic_Empty
     --  then put_line(file,"Corresponding compactification :");
     --       Toric_Compactification_Data.Numeric_Update(-sv);
     --       declare
     --         p : constant Poly_Sys
     --           := Toric_Compactification_Data.Polynomial_System;
     --       begin
     --         put(file,p);
     --         put_line(file,"The face system equals : ");
     --         put(file,Substitute(p'last+1,CMPLX(0.0),p));
     --       end;
     --       Toric_Compactification_Data.Numeric_Clear;
     -- end if;
    end Write_Direction;

  begin
    for i in v'range loop
      Write_Direction(file,v(i).all,errv(i),i);
    end loop;
  end Write_Directions;

end Drivers_for_Path_Directions;
