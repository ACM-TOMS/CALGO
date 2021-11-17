with integer_io;          use integer_io;
with Complex_Numbers_io;  use Complex_Numbers_io;
with Solutions_io;        use Solutions_io;

package body Process_io is

  use Floating_Point_Numbers.double_float_io;

  out_code : output_code;

  procedure Write_path ( ft : in file_type; n : in positive ) is
  begin
    if out_code /= nil
     then put(ft,"*****       path");
          put(ft,n);
          put(ft,"       *****");
          new_line(ft);
    end if;
  end Write_path;

  procedure Write_path ( n : in positive ) is
  begin
    Write_path(Standard_Output,n);
  end Write_path;

  procedure Write_block ( ft : in file_type; n : in positive ) is
  begin
   --if out_code /= nil 
   -- then 
       put(ft,"#####       block");
       put(ft,n);
       put_line(ft,"       #####");
   --end if;
  end Write_block;

  procedure Write_block ( n : in positive ) is
  begin
    Write_Block(Standard_Output,n);
  end Write_block;

  procedure Set_Output_Code ( u : in output_code ) is
  begin
    out_code := u;
  end Set_output_code;

  procedure sWrite ( ft : in file_type; sol : in Solution ) is
  begin
    if (out_code = s) or (out_code = sp) or (out_code = sc) or (out_code = spc)
     then put(ft,sol);
          new_line(ft);
    end if;
  end sWrite;

  procedure sWrite ( sol : in Solution ) is
  begin
    sWrite(Standard_Output,sol);
  end sWrite;

  procedure pWrite ( ft : in file_type; 
                     step : in double_float; t : in double_complex ) is
  begin
    if (out_code = p) or (out_code = sp) or (out_code = pc) or (out_code = spc)
     then put(ft,"step :"); put(ft,step); put(ft,"  ");
          put(ft,"t :"); put(ft,t); new_line(ft);
    end if;
  end pWrite;

  procedure pWrite ( ft : in file_type; step : in double_float;
                     t : in double_complex; sol : in Solution ) is
  begin
    if (out_code = p) or (out_code = sp) or (out_code = pc) or (out_code = spc)
     then put(ft,"step :"); put(ft,step); put(ft,"  ");
          put(ft,"t :"); put(ft,t); new_line(ft);
          if (out_code = sp) or (out_code = spc)
           then put_line(ft,"the predicted solution for t :");
                put_vector(ft,sol);
          end if;
    end if;
  end pWrite;

  procedure pWrite ( step : in double_float; t : in double_complex ) is
  begin
    pWrite(Standard_Output,step,t);
  end pWrite;

  procedure pWrite ( step : in double_float; t : in double_complex;
                     sol : in Solution ) is
  begin
    pWrite(Standard_Output,step,t,sol);
  end pWrite;

  procedure cWrite ( ft : in file_type;
                     normax,normrx,normaf,normrf : in double_float ) is
  begin
    if (out_code = c) or (out_code = pc) or (out_code = sc) or (out_code = spc)
     then put(ft,"correction (a&r):");
          put(ft,normax,3,3,3); put(ft,normrx,3,3,3); put(ft," ");
          put(ft,"residual (a&r):");
          put(ft,normaf,3,3,3); put(ft,normaf,3,3,3); new_line(ft);
    end if;
  end cWrite;

  procedure cWrite ( normax,normrx,normaf,normrf : in double_float ) is
  begin
    cWrite(Standard_Output,normax,normrx,normaf,normrf);
  end cWrite;

  procedure cWrite ( ft : in file_type; 
                     rcond : in double_float; m : in natural ) is
  begin
    if (out_code = c) or (out_code = sc) or (out_code = pc) or (out_code = spc)
     then put(ft,"rcond :"); put(ft,rcond); put(ft,"  ");
          put(ft,"multiplicity : "); put(ft,m,1); new_line(ft);
    end if;
  end cWrite;

  procedure cWrite ( rcond : in double_float; m : in natural ) is
  begin
    cWrite(Standard_Output,rcond,m);
  end cWrite;

  procedure Write_convergence_factor ( factor : in double_float ) is
  begin
    if (out_code = c) or (out_code = sc) or (out_code = pc) or (out_code = spc)
     then put("convergence ratio :"); put(factor); new_line;
    end if;
  end Write_convergence_factor;

  procedure Write_convergence_factor 
                    ( ft : in file_type; factor : in double_float ) is
  begin
    if (out_code = c) or (out_code = sc) or (out_code = pc) or (out_code = spc)
     then put(ft,"convergence ratio :"); put(ft,factor); new_line(ft);
    end if;
  end Write_convergence_factor;

  procedure Write_Statistics ( nstep,nfail,niter,nsyst : in natural ) is
  begin
    Write_Statistics(Standard_Output,nstep,nfail,niter,nsyst);
  end Write_Statistics;

  procedure Write_Statistics ( ft : in file_type;
                               nstep,nfail,niter,nsyst : in natural ) is
  begin
  --if out_code /= nil
  -- then
       put_line(ft,"######################################################");
       put(ft,"number of steps      :"); put(ft,nstep); new_line(ft);
       put(ft,"number of failures   :"); put(ft,nfail);  new_line(ft);
       put(ft,"number of iterations :"); put(ft,niter);   new_line(ft);
       put(ft,"number of systems    :"); put(ft,nsyst); new_line(ft);
  --end if;
  end Write_Statistics;

  procedure Write_Total_Statistics
                        ( tnstep,tnfail,tniter,tnsyst : in natural ) is
  begin
    Write_Total_Statistics(Standard_Output,tnstep,tnfail,tniter,tnsyst);
  end Write_Total_Statistics;

  procedure Write_Total_Statistics 
                ( ft : in file_type;
                  tnstep,tnfail,tniter,tnsyst : in natural ) is
  begin
    put(ft,"total number of steps      :"); put(ft,tnstep); new_line(ft);
    put(ft,"total number of failures   :"); put(ft,tnfail); new_line(ft);
    put(ft,"total number of iterations :"); put(ft,tniter); new_line(ft);
    put(ft,"total number of systems    :"); put(ft,tnsyst); new_line(ft);
  end Write_Total_Statistics;

  procedure sWrite_Solutions ( sols : in Solution_List ) is
  begin
    sWrite_Solutions(Standard_Output,sols);
  end sWrite_Solutions;

  procedure sWrite_Solutions ( ft : in file_type; sols : in Solution_List ) is
  begin
    if out_code /= nil
     then put(ft,sols);
    end if;
  end sWrite_Solutions;

end Process_io;
