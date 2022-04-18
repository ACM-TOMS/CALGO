with integer_io;                         use integer_io;
with File_Operations;                    use File_Operations;
with Timing_Package;                     use Timing_Package;
with Floating_Point_Numbers;             use Floating_Point_Numbers; 
with Float_Vectors_of_Vectors;
with Float_Vectors,Float_Vectors_io;     use Float_Vectors,Float_Vectors_io;
with Integer_Vectors,Integer_Vectors_io; use Integer_Vectors_io;
with Integer_Vectors_of_Vectors;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;      use Complex_Polynomial_Systems_io;

procedure valipoco ( pocofile,resultfile : in file_type ) is

  use Floating_Point_Numbers.double_float_io;

  timer : timing_widget;

-- SCANNING THE INFORMATION FROM FILE :

  procedure Scan_System
                ( file : in file_type; p : out Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   Scans the file for the target system.

    found : boolean := false;
    lp : Link_to_Poly_Sys;

  begin
    Scan_and_Skip(file,"TARGET SYSTEM",found);
    if found
     then get(file,lp); p := lp;
    end if;
  end Scan_System;

  procedure Scan_Dimensions
                ( file : in file_type; n,npaths : out natural ) is

  -- DESCRIPTION :
  --   Scans the input file for the dimension and the number of paths.

    found : boolean := false;
    tmp : integer := 0;

  begin
    Scan_and_Skip(file,"START SOLUTIONS",found);
    if not found
     then n := 0; npaths := 0;
     else get(file,tmp); npaths := tmp;
          get(file,tmp); n := tmp;
    end if;
  end Scan_Dimensions;

  procedure Scan_Data
                ( infile,outfile : in file_type; n,npaths : in natural;
                  nv : out Float_Vectors_of_Vectors.Vector;
                  em : out Integer_Vectors.Vector;
                  ev,rv : out Float_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Scans for the computed directions nv, the errors ev,
  --   the residuals rv and the estimated multiplicities em.

  -- REQUIRED : nv'range = em'range = ev'range = rv'range = 1..npaths.

  -- ON ENTRY :
  --   infile     file with input data from poco;
  --   outfile    file to write results to;
  --   n          dimension of the vectors along the paths;
  --   npaths     number of paths.

  -- ON RETURN :
  --   nv         the computed path directions;
  --   em         estimated multiplicities;
  --   rv         residuals.

    v : Float_Vectors_of_Vectors.Vector(1..npaths);
    m : Integer_Vectors.Vector(1..npaths);
    e,r : Float_Vectors.Vector(1..npaths);
    found : boolean := false;

  begin
    for i in v'range loop                                -- scan for normals
      Scan_and_Skip(infile,"Computed direction",found);
      exit when not found;
      v(i) := new Float_Vectors.Vector(1..n);
      get(infile,v(i).all);
      Scan(infile,"error :",found);
      exit when not found;
      get(infile,e(i));
    end loop;
    for i in m'range loop                  -- scan for normals and residuals
      Scan(infile,"m :",found);
      exit when not found;
      get(infile,m(i));
      Scan(infile,"res :",found);
      exit when not found;
      get(infile,r(i));
    end loop;
    put_line(outfile,
             "COMPUTED DIRECTIONS, ESTIMATED MULTIPLICITY AND RESIDUALS : ");
    for i in v'range loop
      put(outfile,i,1); put(outfile," :"); put(outfile,v(i),3,3,3);
      put(outfile," err : "); put(outfile,e(i),3,3,3);
      put(outfile," m : "); put(outfile,m(i),1);
      put(outfile," res : "); put(outfile,r(i),2,3,3);
      new_line(outfile);
    end loop;
    nv := v; em := m; ev := e; rv := r;
  end Scan_Data;

-- PROCESSING THE DATA :

  function Lattice_Error ( x : Float_Vectors.Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the sum of the differences of the components in x with
  --   the closest corresponding integer.

    res : double_float := 0.0;

  begin
    for i in x'range loop
      res := res + abs(x(i) - double_float(integer(x(i))));
    end loop;
    return res;
  end Lattice_Error;

  function Lattice_Errors ( v : Float_Vectors_of_Vectors.Vector )
                          return Float_Vectors.Vector is

    res : Float_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Lattice_Error(v(i).all);
    end loop;
    return res;
  end Lattice_Errors;

  function Maximum ( v : Float_Vectors.Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the largest component in absolute value.

    max : double_float := abs(v(v'first));

  begin
    for i in v'first+1..v'last loop
      if abs(v(i)) > max
       then max := abs(v(i));
      end if;
    end loop;
    return max;
  end Maximum;

  procedure Scale ( v : in out Float_Vectors_of_Vectors.Vector;
                    e : in out Float_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Divides every vector by its largest element in absolute value,
  --   with as well the corresponding errors.

    use Float_Vectors;
    tol : double_float := 10.0**(-8);
    max : double_float;

  begin
    for i in v'range loop
      max := Maximum(v(i).all);
      if max > tol
       then v(i).all := (1.0/max)*v(i).all;
            e(i) := e(i)/max;
      end if;
    end loop;
  end Scale;

-- COMPUTE FREQUENCY TABLE :

  function Equal ( v1,v2 : Float_Vectors.Vector; tol : double_float )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if abs(v1(i)-v2(i)) <= tol, for i in v1'range=v2'range.

  begin
    for i in v1'range loop
      if abs(v1(i)-v2(i)) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  procedure Update_Frequency 
              ( vpath : in out Integer_Vectors_of_Vectors.Vector;
                i,pathno : in natural ) is

  -- DESCRIPTION :
  --   A path, with pathno, has as path direction i.

    freq : Integer_Vectors.Vector(vpath(i)'first..vpath(i)'last+1);

  begin
    freq(vpath(i)'range) := vpath(i).all;
    freq(freq'last) := pathno;
    Integer_Vectors.Clear(vpath(i));
    vpath(i) := new Integer_Vectors.Vector'(freq);
  end Update_Frequency;

  procedure Update_Frequency_Table
              ( freqv : in out Float_Vectors_of_Vectors.Vector;
                vpath : in out Integer_Vectors_of_Vectors.Vector;
                cnt : in out natural; tol : in double_float;
                v : in Float_Vectors.Vector; pathno : in natural ) is

  -- DESCRIPTION :
  --   Scans the current frequency table for the vector v.

    done : boolean := false;

  begin
    for i in 1..cnt loop
      if Equal(v,freqv(i).all,tol)
       then Update_Frequency(vpath,i,pathno); done := true;
      end if;
      exit when done;
    end loop;
    if not done
     then cnt := cnt + 1;
          freqv(cnt) := new Float_Vectors.Vector'(v);
          vpath(cnt) := new Integer_Vectors.Vector'(1..1 => pathno);
    end if;
  end Update_Frequency_Table;

  procedure Frequency_Table
              ( v : in Float_Vectors_of_Vectors.Vector;
                e,r : in Float_Vectors.Vector; tol : in double_float;
                freqv : out Float_Vectors_of_Vectors.Link_to_Vector;
                vpath : out Integer_Vectors_of_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Computes the frequency table of path directions.
  --   Only those directions for which the error is smaller than tol
  --   and the corresponding residual is higher than tol are considered.

  -- ON ENTRY :
  --   v        scaled computed path directions;
  --   e        errors on the computed path directions;
  --   r        residuals of the solutions at the end of the paths;
  --   tol      tolerance for precision.

  -- ON RETURN :
  --   freqv    different path directions;
  --   vpath    for each direction, vector of paths with same direction.

    cnt : natural := 0;
    freqdirs : Float_Vectors_of_Vectors.Vector(v'range);
    freqpath : Integer_Vectors_of_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      if r(i) > tol and then e(i) < tol
       then Update_Frequency_Table(freqdirs,freqpath,cnt,tol,v(i).all,i);
      end if;
    end loop;
    freqv := new Float_Vectors_of_Vectors.Vector'(freqdirs(1..cnt));
    vpath := new Integer_Vectors_of_Vectors.Vector'(freqpath(1..cnt));
  end Frequency_Table;

  function Average_Error ( pathdir : Integer_Vectors.Vector;
                           e : Float_Vectors.Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the average error over all paths in pathdir.

    cnt : constant natural := pathdir'length;
    res : double_float := 0.0;

  begin
    for i in pathdir'range loop
      res := res + e(pathdir(i));
    end loop;
    res := res/double_float(cnt);
    return res;
  end Average_Error;

  procedure Write_Frequency_Table 
               ( file : in file_type;
                 freqv : in Float_Vectors_of_Vectors.Vector;
                 vpath : in Integer_Vectors_of_Vectors.Vector;
                 m : in Integer_Vectors.Vector; e : in Float_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the frequency table on file, together with estimated m and error.

  begin
    put_line(file,"FREQUENCY TABLE OF PATH DIRECTIONS :");
    for i in freqv'range loop
      put(file,"v("); put(file,i,1); put(file,") : ");
      put(file,freqv(i),3,3,3); 
      put(file,"  m : "); put(file,m(vpath(i)(1)),1);
      put(file," err : "); put(file,e(vpath(i)(1)),3,3,3);
      put(file," avgerr : "); put(file,Average_Error(vpath(i).all,e),3,3,3);
      put(file," laterr : "); put(file,Lattice_Error(freqv(i).all),3,3,3);
      new_line(file);
      put(file," with "); put(file,vpath(i)'last,1);
      put(file," paths : "); put(file,vpath(i)); new_line(file);
    end loop;
  end Write_Frequency_Table;

-- COMPUTING FACES OF POLYNOMIAL SYSTEMS :

  function Product ( d : Degrees; v : Float_Vectors.Vector )
                   return double_float is

  -- DESCRIPTION : Returns <d,v>.
 
    res : double_float := 0.0;

  begin
    for i in d'range loop
      res := res + double_float(d(i))*v(i);
    end loop;
    return res;
  end Product;

  function Minimal_Support ( p : Poly; v : Float_Vectors.Vector )
                           return double_float is

    min : double_float := 0.0;
    first : boolean := true;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is

      ip : double_float := Product(t.dg,v);

    begin
      if first 
       then min := ip; first := false;
       elsif ip < min
           then min := ip;
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return min;
  end Minimal_Support;

  function Face ( p : Poly; tol,ip : double_float; v : Float_Vectors.Vector )
                return Poly is

  -- DESCRIPTION :
  --   Returns those terms t for which abs(<t.dg,v> - ip) <= tol.

    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is
    begin
      if abs(Product(t.dg,v) - ip) <= tol
       then Plus_Term(res,t);
      end if;
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Face;

  function Face ( p : Poly; tol : double_float; v : Float_Vectors.Vector ) 
                return Poly is

  -- DESCRIPTION :
  --   Returns the face of the polynomial with outer normal v.

    ip : double_float := Minimal_Support(p,v);

  begin
    return Face(p,tol,ip,v);
  end Face;

  function Face ( p : Poly_Sys; tol : double_float; v : Float_Vectors.Vector )
                return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the face of the system with outer normal v.

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Face(p(i),tol,v);
    end loop;
    return res;
  end Face;

  procedure Face_Systems
               ( file : in file_type; p : in Poly_Sys; tol : in double_float;
                 v : in Float_Vectors_of_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the faces and the number of paths that diverged to it.

  begin
    put_line(file,"FACE SYSTEMS :");
    for i in v'range loop
      put(file,"path direction ");
      put(file,i,1); put_line(file," :");
      put(file,Face(p,tol,v(i).all));
    end loop;
  end Face_Systems;

-- MAIN PROCEDURE :

  procedure Main ( infile,outfile : in file_type ) is

    tol : constant double_float := 10.0**(-2);
    lp : Link_to_Poly_Sys;
    n,npaths : natural;

  begin
    Scan_System(infile,lp);
    put_line(outfile,"TARGET SYSTEM : "); put(outfile,lp.all);
    Scan_Dimensions(infile,n,npaths);
    put(outfile," n : "); put(outfile,n,1);
    put(outfile," and #paths : "); put(outfile,npaths,1); new_line(outfile);
    declare
      v : Float_Vectors_of_Vectors.Vector(1..npaths);
      m : Integer_Vectors.Vector(1..npaths);
      e,r,le : Float_Vectors.Vector(1..npaths);
      freqv : Float_Vectors_of_Vectors.Link_to_Vector;
      vpath : Integer_Vectors_of_Vectors.Link_to_Vector;
    begin
      Scan_Data(infile,outfile,n,npaths,v,m,e,r);
      le := Lattice_Errors(v);
     -- Scale(v,e); -- avoid because, then fractions may arise
      Frequency_Table(v,e,r,tol,freqv,vpath);
      Write_Frequency_Table(outfile,freqv.all,vpath.all,m,e);
      Face_Systems(outfile,lp.all,tol,freqv.all);
    end;
  end Main;

begin
  tstart(timer);
  Main(pocofile,resultfile);
  tstop(timer);
  new_line(resultfile);
  print_times(resultfile,timer,"Validation stage of polyhedral end game");
  new_line(resultfile);
end valipoco;
