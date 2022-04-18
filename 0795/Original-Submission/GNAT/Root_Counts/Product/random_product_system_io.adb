with integer_io;
with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers,Complex_Vectors;
with Natural_Vectors;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Symbol_Table,Symbol_Table_io;
with Complex_Multivariate_Polynomials_io;
 use Complex_Multivariate_Polynomials_io;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;      use Complex_Polynomial_Systems_io;
with Random_Product_System;

package body Random_Product_System_io is

  use Floating_Point_Numbers.double_float_io;

  procedure get ( n : out natural ) is
  begin
    Random_Product_System_io.get(Standard_Input,n);
  end get;

  procedure get ( file : in file_type; n : out natural ) is

    nn : natural;

  begin
    integer_io.get(file,nn); n := nn;
    Random_Product_System.Init(nn);
    declare
      h : Vector(0..nn);
      pp : Poly;
      d : degrees := new Natural_Vectors.Vector'(1..nn => 0);
      stop : boolean;
    begin
      get(file,nn,pp); Clear(pp);
      for i in 1..nn loop
        stop := false;
        while not stop loop
          get(file,nn,pp);
          stop := (pp = Null_Poly);
          exit when stop;
          h(0) := Coeff(pp,d);
          for j in 1..nn loop
            d(j) := 1;
            h(j) := Coeff(pp,d);
            d(j) := 0;
          end loop;
          Random_Product_System.Add_Hyperplane(i,h);
        end loop;
      end loop;
    end;
  end get;

  procedure put ( n,fore,after,exp : in natural ) is
  begin
    Random_Product_System_io.put(Standard_Output,n,fore,after,exp);
  end put;

  procedure put ( file : in file_type; n,fore,after,exp : in natural ) is

    h : Vector(0..n);

    procedure Write_Number ( file : in file_type; x : in double_complex ) is
    begin
      if IMAG_PART(x) + 1.0 = 1.0
       then put(file,REAL_PART(x),fore,after,exp);
       else put(file,'(');
            put(file,REAL_PART(x),fore,after,exp);
            put(file,'+');
            put(file,IMAG_PART(x),fore,after,exp);
            put(file,')');
      end if;
    end Write_Number;

  begin
    for i in 1..n loop
      put(file,"The hyperplanes for the "); integer_io.put(file,i,1); 
      put_line(file,"th equation :");
      for j in 1..Random_Product_System.Number_Of_Hyperplanes(i) loop
        h := Random_Product_System.Get_Hyperplane(i,j);
        put(file,' ');
        for k in 1..n loop
  	  Write_Number(file,h(k));
	  put(file,'*');
	  Symbol_Table_io.put(file,Symbol_Table.Get(k));
	  put(file," + ");
        end loop;
        Write_Number(file,h(0));
        new_line(file);
      end loop;
    end loop;
  end put;

end Random_Product_System_io;
