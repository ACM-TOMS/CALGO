with unchecked_deallocation;
with Lists;
with Random_Product_System;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers;
with Random_Number_Generators;

package body Templates is

 -- DATA STRUCTURES :

  package List_of_Vectors is new Lists(Link_to_Vector);
  type Equation_List is new List_of_Vectors.List;
  
  type Equation is record
    first,last : Equation_List;
  end record;

  type Equations is array(positive range <>) of Equation;
  type Link_To_Equations is access Equations;
  procedure free is new unchecked_deallocation (Equations,Link_To_Equations);

 -- INTERNAL DATA :
 
  rps : Link_To_Equations;

--------------------
--  CONSTRUCTORS  --
--------------------

  procedure Create ( n : in natural ) is
  begin
    rps := new Equations(1..n);
  end Create;

  procedure Add_Hyperplane ( i : in natural; h : in Vector ) is

    eqi : Equation renames rps(i);
    lh : Link_To_Vector := new Vector'(h);

  begin
    if Is_Null(eqi.first)
     then Construct(lh,eqi.first);
          eqi.last := eqi.first;
     else declare 
            temp : Equation_List;
          begin
            Construct(lh,temp);
            Swap_Tail(eqi.last,temp);
            eqi.last := Tail_Of(eqi.last);
          end;
    end if;
  end Add_Hyperplane;

  procedure Change_Hyperplane ( i,j : in natural; h : in Vector ) is
  begin
    if rps = null
     then return;
     else declare
            eqi : Equation_List := rps(i).first;
	    lv : Link_To_Vector;
            count : natural := 1;
          begin
            while not Is_Null(eqi) loop
              if count = j
               then lv := Head_Of(eqi);
		    for k in h'range loop
		      lv(k) := h(k);
                    end loop;
                    return;
               else count := count + 1;
                    eqi := Tail_Of(eqi);
              end if;
            end loop;
          end;
    end if;
  end Change_Hyperplane;

-----------------
--  SELECTORS  --
-----------------

  function Number_of_Hyperplanes ( i : natural ) return natural is
  begin
    if rps = null
     then return 0;
     else return Length_Of(rps(i).first);
    end if;
  end Number_of_Hyperplanes;

  procedure Get_Hyperplane ( i,j : in natural; h : out Vector ) is
  begin
    h := (h'range => 0);
    if rps = null
     then return;
     else declare
            eqi : Equation_List := rps(i).first;
            count : natural := 1;
          begin
            while not Is_Null(eqi) loop
              if count = j
               then h := Head_Of(eqi).all;
                    return;
               else count := count + 1;
                    eqi := Tail_Of(eqi);
              end if;
            end loop;
          end;
    end if;
  end Get_Hyperplane;

  procedure Polynomial_System ( n,nbfree : in natural ) is

    rndms : Complex_Vectors.Vector(0..nbfree);

  begin
   -- GENERATE THE FREE COEFFICIENTS :
    rndms(0) := CMPLX(0.0);
    for i in rndms'first+1..rndms'last  loop
     -- rndms(i) := CMPLX(Random_Number_Generators.Random);
      rndms(i) := Random_Number_Generators.Random1;
       -- random complex number with modulus one
    end loop;
   -- BUILD THE RANDOM PRODUCT SYSTEM :
    Random_Product_System.Init(n);
    for i in 1..n loop
      for j in 1..Number_of_Hyperplanes(i) loop
        declare
	  ih : Natural_Vectors.Vector(0..n);
	  h : Complex_Vectors.Vector(0..n);
        begin
	  Get_Hyperplane(i,j,ih);
	  for k in h'range loop
	    h(k) := rndms(ih(k));
          end loop;
          Random_Product_System.Add_Hyperplane(i,h);
        end;
      end loop;
    end loop;
  end Polynomial_System;

  function Verify ( n : natural; lp : List ) return natural is

    temp : List := lp;
    stop : boolean := false;        
    matrix : array (1..n,1..n) of natural;
    nb : natural;

    function Degenerate return boolean is
      degen : boolean;
      first : natural;
    begin
      for i in 1..n loop
	first := matrix(i,1);
	degen := true;
	for j in 2..n loop
	  if matrix(i,j) /= first
	   then degen := false;
          end if;
	  exit when not degen;
        end loop;
	if degen
	 then return true;
        end if;
      end loop;
      for j in 1..n loop
	first := matrix(1,j);
	degen := true;
	for i in 2..n loop
	  if matrix(i,j) /= first
	   then degen := false;
          end if;
	  exit when not degen;
        end loop;
	if degen
	 then return true;
        end if;
      end loop;
      return false;
    end Degenerate;

    procedure PVerify ( i,n : in natural; sum : in out natural ) is
    begin
      if i > n
       then if Is_Null(temp)
	     then sum := sum + 1;
		  stop := true;
	     elsif Degenerate
		 then stop := true;
		 else temp := Tail_Of(temp);
		      sum := sum + 1;
            end if;
       else declare
              eqi : Equation_List := rps(i).first;
              h : Vector(0..n);
              count : natural := 0;
             begin
               while not Is_Null(eqi) loop
                 count := count + 1;
    	         if count = Head_Of(temp)(i)
                  then h := Head_Of(eqi).all;
                       for j in 1..n loop
                         matrix(i,j) := h(j);
                       end loop;
                       PVerify(i+1,n,sum);
                 end if;
	       	 exit when stop;
                 eqi := Tail_Of(eqi);
               end loop;
             end;
       end if;
     end PVerify;

  begin
    nb := 0;
    if not Is_Null(temp)
     then PVerify(1,n,nb);
    end if;
    return nb;
  end Verify;

------------------
--  DESTRUCTOR  --
------------------

  procedure Clear ( eql : in out Equation_List ) is

    temp : Equation_List := eql;
    lv : Link_To_Vector;

  begin
    while not Is_Null(temp) loop
      lv := Head_Of(temp);
      Clear(lv);
      temp := Tail_of(temp);
    end loop;
    List_Of_Vectors.Clear(List_Of_Vectors.List(eql));
  end Clear;

  procedure Clear ( eq : in out Equation ) is
  begin
    Clear(eq.first);
    -- eq.last is just a pointer to the last element of eq.first;
    -- if eq.first disappears, then also eq.last does
  end Clear;

  procedure Clear ( eqs : in out Equations ) is
  begin
    for i in eqs'range loop
      Clear(eqs(i));
    end loop;
  end Clear;

  procedure Clear is
  begin
    if rps /= null
     then for i in rps'range loop
            Clear(rps(i));
          end loop;
          free(rps);
    end if;
  end Clear;

end Templates;
