with unchecked_deallocation;
with Complex_Numbers;                  use Complex_Numbers;
with Greatest_Common_Divisor;          use Greatest_Common_Divisor;
with Integer_Linear_System_Solvers;    use Integer_Linear_System_Solvers;

package body Transformations is

  type Transfo_tp is new matrix;
  procedure free is new unchecked_deallocation (Transfo_tp,Transfo);

-- CONSTRUCTORS :

  function Id ( n : natural ) return Transfo is
    t : Transfo;
  begin
    t := new Transfo_tp(1..n,1..n);
    for i in 1..n loop
      for j in 1..n loop
	t(i,j) := 0;
      end loop;
      t(i,i) := 1;
    end loop;
    return t;
  end Id;

  function Create ( v : Vector; i : integer ) return Transfo is
    t : Transfo;
    j : integer;
  begin
    t := Id(v'last-v'first+1);
    if v(i) = 0
     then return t;
     else j := 0;
	  for k in v'range loop
	    if (v(k) /= 0) and (k /= i)
             then j := k;
            end if;
	    exit when j /= 0;
          end loop;
	  if j = 0
	   then return t;
	   else declare
		  a,b,k,l,d : integer;
                begin
		  a := v(i); b := v(j);
		  gcd(a,b,k,l,d);
		  a := a/d;  b := b/d;
		  t(i,i) := k;  t(i,j) := l;
		  t(j,i) := -b; t(j,j) := a;
                end;
		return t;
          end if;
    end if;
  end Create;

  function Create ( v : Integer_Vectors_of_Vectors.Vector ) return Transfo is
    t : Transfo;
  begin
    t := new Transfo_tp(v'range,v'range);
    for i in v'range loop
      for j in v(i)'range loop
	t(j,i) := v(i)(j);
      end loop;
    end loop;
    return t;
  end Create;

  function Create ( m : Matrix ) return Transfo is
    t : Transfo;
  begin
    t := new Transfo_tp(m'range(1),m'range(2));
    t.all := transfo_tp(m);
    return t;
  end Create;

  function Rotate ( v : Vector; i : integer ) return Transfo is
    t : Transfo;
    j : integer;
  begin
    t := Id(v'last-v'first+1);
    if v(i) = 0
     then return t;
     else declare
	    w : Vector(v'range) := v;
	    acc : Transfo := new Transfo_tp'(t.all);
	  begin
	    loop
  	      j := 0;
	      for k in w'range loop
	        if (w(k) /= 0) and (k /= i)
                 then j := k;
                end if;
	        exit when j /= 0;
              end loop;
	      if j /= 0
	       then declare
		      a,b,k,l,d : integer;
                    begin
		      a := w(i); b := w(j);
		      gcd(a,b,k,l,d);
		      a := a/d;  b := b/d;
		      t(i,i) := k;  t(i,j) := l;
		      t(j,i) := -b; t(j,j) := a;
                    end;
		    Apply(t,w);
		    Mult2(t,acc);
		    t(i,i) := 1; t(j,j) := 1;
		    t(i,j) := 0; t(j,i) := 0;
              end if;
	      exit when j = 0;
            end loop;
	    Clear(t);
	    return acc;
          end;
    end if;
  end Rotate;

  function Rotate ( v : Link_to_Vector; i : integer ) return Transfo is
  begin
    if v = null
     then return null;
     else return Rotate(v.all,i);
    end if;
  end Rotate;

  procedure Transpose ( t : in out Transfo ) is

   -- DESCRIPTION :
   --   Transposes the matrix of t.

    temp : integer;

  begin
    for i in t'range(1) loop
      for j in t'first(2)..(i-1) loop
	temp := t(i,j);
	t(i,j) := t(j,i);
	t(j,i) := temp;
      end loop;
    end loop;
  end Transpose;

  function Build_Transfo ( v : Vector; i : integer ) return Transfo is
    t : Transfo;
    ind : integer;
    zeroes : boolean;
  begin
    t := Id(v'last-v'first+1);
    if v(i) /= 0
     then ind := i;
     else ind := v'last + 1;
          for k in v'range loop
            if v(k) /= 0
             then ind := k;
	          exit;
            end if;
          end loop;
    end if;
    if ind = v'last + 1
     then return t;
     else declare
	    w : Vector(v'range) := v;
	    acc : Transfo := new Transfo_tp'(t.all);
	    j1,j2 : integer;
	  begin
	    if v(i) = 0
	     then acc(i,i) := 0; acc(i,ind) := 1;
		  acc(ind,i) := 1; acc(ind,ind) := 0;
		  w(i) := w(ind); w(ind) := 0;
            end if;
            for k in v'first..i loop
	      if w(k) /= 0
	       then j1 := k;
		    exit;
              end if;
            end loop;
	    if j1 /= i
             then zeroes := false;
		  loop
		    for k in (j1+1)..i loop
		      if w(k) /= 0
		       then j2 := k;
		            exit;
                      end if;
                    end loop;
	            declare
	              a,b,k,l,d : integer;
                    begin
	              a := w(j1); b := w(j2);
	              gcd(a,b,k,l,d);
	              a := a/d;  b := b/d;
	              t(j1,j1) := l;  t(j1,j2) := -k;
	              t(j2,j1) := a; t(j2,j2) := b;
	              w(j2) := d;
                    end;
	            Mult2(t,acc);
  	            t(j1,j1) := 1; t(j2,j2) := 1;
	            t(j1,j2) := 0; t(j2,j1) := 0;
		    if j2 < i
		     then j1 := j2;
		     else exit;
                    end if;
                  end loop;
             else zeroes := true;
            end if;
            for k in reverse i..v'last loop
	      if w(k) /= 0
	       then j2 := k;
		    exit;
              end if;
            end loop;
	    if j2 /= i
             then loop
		    for k in reverse i..(j2-1) loop
		      if w(k) /= 0
		       then j1 := k;
		            exit;
                      end if;
                    end loop;
	            declare
	              a,b,k,l,d : integer;
                    begin
	              a := w(j1); b := w(j2);
	              gcd(a,b,k,l,d);
	              a := a/d;  b := b/d;
	              t(j1,j1) := a;  t(j1,j2) := b;
	              t(j2,j1) := -l; t(j2,j2) := k;
	              w(j1) := d;
                    end;
	            Mult2(t,acc);
  	            t(j1,j1) := 1; t(j2,j2) := 1;
	            t(j1,j2) := 0; t(j2,j1) := 0;
		    if j1 > i
		     then j2 := j1;
		     else exit;
                    end if;
                  end loop;
             elsif zeroes
		 then if w(i) < 0
		       then t(i,i) := -1;
			    Mult2(t,acc);
                      end if;
	    end if;
	    Clear(t); --otherwise segmentation fault...
	    return acc;
          end;
    end if;
  end Build_Transfo;

  function Build_Transfo ( v : Link_to_Vector; i : integer ) return Transfo is
  begin
    if v = null
     then return null;
     else return Build_Transfo(v.all,i);
    end if;
  end Build_Transfo;

-- SELECTOR :

  function Dimension ( t : Transfo ) return natural is
  begin
    return (t'last(1) - t'first(1) + 1);
  end Dimension;

  function Sign ( t : Transfo ) return integer is
  begin
    if t /= null
     then declare
            a  : matrix(t'range(1),t'range(2)) := matrix(t.all);
          begin
            return Det(a);
          end;
     else return 0;
    end if;
  end Sign;

-- OPERATIONS :

  function Transpose ( t : Transfo ) return Transfo is

    res : Transfo := new Transfo_tp(t'range(2),t'range(1));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := t(j,i);
      end loop;
    end loop;
    return res;
  end Transpose;

  function Invert ( t : Transfo ) return Transfo is

    n : constant natural := Dimension(t);
    triang : matrix(t'range(1),t'range(2)) := matrix(t.all); 
    l : matrix(t'range(1),t'range(1));
    res : Transfo := new Transfo_tp(t'range(1),t'range(2));
    x,b : Vector(1..n) := (1..n => 0);
    fail : boolean;

  begin
    Upper_Triangulate(l,triang);
    for j in t'range(2) loop
      for i in l'range(1) loop  -- image of jth basis vector
	b(i) := l(i,j);
      end loop;
      Solve1(triang,x,b,fail);
      for i in t'range(1) loop
	res(i,j) := x(i);
      end loop;
    end loop;
    return res;
  end Invert;

  function "*" ( t1,t2 : Transfo ) return Transfo is
    r : Transfo;
  begin
    r := new Transfo_tp(t1'range(1),t2'range(2));
    r.all := t1.all*t2.all;
    return r;
  end "*";

  procedure Mult1 ( t1 : in out Transfo; t2 : in Transfo ) is
  begin
    Mult1(t1.all,t2.all);
  end Mult1;

  procedure Mult2 ( t1 : in Transfo; t2 : in out Transfo ) is
  begin
    Mult2(t1.all,t2.all);
  end Mult2;

  function "*" ( t : Transfo; v : Vector ) return Vector is
    r : Vector(t'range(1));
  begin
    r := t.all*v;
    return r;
  end "*";

  function "*" ( t : Transfo; v : Link_to_Vector ) return Link_to_Vector is
  begin
    if v = null
     then return v;
     else declare
            res : Link_to_Vector := new Vector(t'range(1));
          begin
            res.all := t.all*v.all;
            return res;
          end;
    end if;
  end "*";

  procedure Apply ( t : in Transfo; v : in out Vector ) is
  begin
    v := t.all*v;
  end Apply;

  procedure Apply ( t : in Transfo; v : in out Link_to_Vector ) is
  begin
    if v /= null
     then Apply(t,v.all);
    end if;
  end Apply;

  function "*" ( t : Transfo; v : Complex_Vectors.Vector )
               return Complex_Vectors.Vector is

    res : Complex_Vectors.Vector(v'range);

  begin
    for j in t'range(2) loop
      res(j) := CMPLX(1.0);
      for i in t'range(1) loop
        res(j) := res(j)*v(i)**t(i,j);
      end loop;
    end loop;
    return res;
  end "*";

  procedure Apply ( t : in Transfo; v : in out Complex_Vectors.Vector ) is
  begin
    v := t*v;
  end Apply;

-- DESTRUCTOR :

  procedure Clear ( t : in out Transfo ) is
  begin
    free(t);
  end Clear;

end Transformations;
