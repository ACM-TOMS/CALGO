package body Greatest_Common_Divisor is

  function pos_gcd ( a,b : integer ) return integer is
    r : integer;
  begin
    if b = 0
     then return a;
     else r := a mod b;
          if r = 0
           then return b;
           else return pos_gcd(b,r);
          end if;
    end if;
  end pos_gcd;

  function gcd ( a,b : integer ) return integer is
  begin
    if a < 0
     then if b < 0
	   then return pos_gcd(-a,-b);
	   else return pos_gcd(-a,b);
          end if;
     elsif b < 0
	 then return pos_gcd(a,-b);
	 else return pos_gcd(a,b);
    end if;
  end gcd;

  function lcm ( a,b : integer ) return integer is
    g : integer := gcd(a,b);
  begin
    if g = 0
     then return 0;
     else return (a/g)*b;
    end if;
  end lcm;

  procedure gcd ( a,b : in integer; k,l,d : out integer ) is

    kk,ll : integer;

    procedure pos_gcd ( a,b : in integer; k,l,d : out integer ) is

      r,aa,bb,k1,l1,k0,l0,h,tmp : integer;

     -- REQUIRED :
     --   a >= b and a > 0 and b > 0

    begin
      if a mod b = 0
       then k := 0; l := 1; d := b;
       else k0 := 0; l0 := 1;
            r := a mod b;
            aa := b;
            bb := r;
            k1 := 1; l1 := - a / b;
            loop
              r := aa mod bb;
              exit when r = 0;
              h := aa / bb;
  	      tmp := k1;
  	      k1 := k0 - h*k1; 
	      k0 := tmp;
	      tmp := l1;
	      l1 := l0 - h*l1;
	      l0 := tmp;
	      aa := bb;
	      bb := r;
            end loop;
	    k := k1; l := l1; d := bb;
      end if;
    end pos_gcd;

  begin
    if a = 0
     then if b < 0
           then d := -b; k := 0; l := -1;
           else d := b; k := 0; l := 1;
          end if;
	  return;
     elsif b = 0
         then if a < 0
               then d := -a; k := -1; l := 0;
               else d := a; k := 1; l := 0;
              end if;
	      return;
    end if;
    if a < b
     then if b < 0
           then pos_gcd(-b,-a,ll,kk,d);
                k := -kk; l := -ll;
  	   elsif a < 0
	       then pos_gcd(b,-a,l,kk,d);
	            k := -kk;
	       else pos_gcd(b,a,l,k,d);
          end if;
     else -- a >= b
          if a < 0
           then pos_gcd(-a,-b,kk,ll,d);
                k := -kk; l := -ll;
           elsif b < 0
               then pos_gcd(a,-b,k,ll,d);
          	    l := -ll;
  	       else pos_gcd(a,b,k,l,d);
          end if;
    end if;
  end gcd;

end Greatest_Common_Divisor;
