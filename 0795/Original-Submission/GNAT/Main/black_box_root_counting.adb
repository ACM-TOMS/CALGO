with integer_io;                         use integer_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Integer_Vectors;                    use Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;      use Complex_Polynomial_Systems_io;
with Solutions_io;                       use Solutions_io;
with Total_Degree_Start_Systems;         use Total_Degree_Start_Systems;
with Partitions_of_Sets_Of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_Of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with m_Homogeneous_Bezout_Numbers;       use m_Homogeneous_Bezout_Numbers;
with m_Homogeneous_Start_Systems;        use m_Homogeneous_Start_Systems;
with Set_Structure,Set_Structure_io;
with Degree_Sets_Tables;                 use Degree_Sets_Tables;
with Random_Product_System;
with Construct_Random_Product_Start_System;
 use Construct_Random_Product_Start_System;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Black_Box_Mixed_Volume_Computations;
 use Black_Box_Mixed_Volume_Computations;

procedure Black_Box_Root_Counting
               ( file : in out file_type;
                 p : in Poly_Sys; rc : out natural;
                 q : out Poly_Sys; qsols : out Solution_List;
                 rocotime,hocotime : out duration ) is

  function Set_Structure_Bound ( p : Poly_Sys ) return natural is
  begin
    Construct_Random_Product_Start_System.Build_Set_Structure(p);
    return Permanent(Degree_Sets_Tables.Create);
  end Set_Structure_Bound;

  procedure Count_Roots 
               ( file : in file_type; p : in Poly_Sys;
                 tode,mhbz,setb,mivo : out natural;
                 zz : out Partition; nz : out natural;
                 lifsup : out Link_to_Array_of_Lists;
                 mix : out Link_to_Vector; mixsub : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Computes four different root counts for the system p.

  -- ON ENTRY :
  --   file      output file;
  --   p         a polynomial system.

  -- ON RETURN :
  --   tode      total degree;
  --   mhbz      m-homogeneous Bezout number;
  --   setb      bound based on set structure;
  --   mivo      mixed volume;
  --   zz        partition used to compute mhbz;
  --   nz        number of sets in partition zz;
  --   lifsup    lifted supports of the system;
  --   mix       type of mixture;
  --   mixsub    mixed subdivision used to compute mivo.

    timer : timing_widget;
    n : constant natural := p'length;
    d,bz,m,bs,mv : natural;
    z : Partition(1..n);

  begin
    tstart(timer);
    d := Total_Degree(p);
    PB(p,bz,m,z);
    bs := Set_Structure_Bound(p);
    Black_Box_Mixed_Volume_Computation(p,mix,lifsup,mixsub,mv);
    tstop(timer);
    new_line(file);
    put_line(file,"ROOT COUNTS :");
    new_line(file);
    put(file,"total degree : ");
      put(file,d,1); new_line(file);
    put(file,m,1); put(file,"-homogeneous Bezout number : ");
      put(file,bz,1); new_line(file);
      put(file,"  with partition : "); put(file,z(1..m)); new_line(file);
    put(file,"generalized Bezout number : ");
      put(file,bs,1); new_line(file);
      put_line(file,"  based on the set structure :");
      Set_Structure_io.put(file);
    put(file,"mixed volume : "); put(file,mv,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Root Counting");
    tode := d; mhbz := bz; setb := bs; mivo := mv;
    zz := z; nz := m;
    rocotime := Elapsed_User_Time(timer);
  end Count_Roots;

  function Is_Minimum ( a,b,c,x : natural ) return boolean is

  -- DESCRIPTION : 
  --   Returns true if x <= a, x <= b, and x <= c.

  begin
    return (x <= a) and (x <= b) and (x <= c);
  end Is_Minimum;

  procedure Construct_Start_System
                 ( file : in file_type; p : in Poly_Sys;
                   d,bz,bs,mv : in natural; z : in Partition;
                   mix : in Vector; lifted : in Array_of_Lists;
                   mixsub : in Mixed_Subdivision; roco : out natural;
                   qq : in out Poly_Sys; qqsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Constructs a start system for the minimal root count and least
  --   amount of work.

  -- ON ENTRY :
  --   file        output file;
  --   p           polynomial system;
  --   d           total degree;
  --   bz          m-homogeneous Bezout number;
  --   bs          Bezout number based on set structure;
  --   mv          mixed volume;
  --   z           partition that corresponds with bz;
  --   mix         type of mixture of the supports;
  --   mixsub      mixed subdivision used to computed mv.

  -- ON RETURN :
  --   roco        minimum(d,bz,bs,mv);
  --   qq          start system;
  --   qqsols      solutions of qq.

    timer : timing_widget;
    n : constant natural := p'length;
    nl : natural;

  begin
    new_line(file);
    tstart(timer);
    if Is_Minimum(bz,bs,mv,d)
     then put_line(file,"START SYSTEM BASED ON TOTAL DEGREE :");
          roco := d;
          Start_System(p,qq,qqsols);
     elsif Is_Minimum(d,bs,mv,bz)
         then put(file,z'length,1);
              put_line(file,"-HOMOGENEOUS START SYSTEM :");
              roco := bz;
              m_Homogeneous_Start_System(p,z,qq,qqsols);
         elsif Is_Minimum(d,bz,mv,bs)
            then put_line(file,"LINEAR-PRODUCT START SYSTEM : ");
                 roco := bs;
                 Random_Product_System.Init(n);
                 Build_Random_Product_System(n);
                 qq := Random_Product_System.Polynomial_System;
                 Random_Product_System.Solve(qqsols,nl);
                 Set_Structure.Clear;
                 Random_Product_System.Clear;
            else put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
                 Black_Box_Polyhedral_Continuation
                   (p,mix,lifted,mixsub,qq,qqsols);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,qq);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line(file);
    put(file,Length_Of(qqsols),Head_Of(qqsols).n,qqsols);
    new_line(file);
    print_times(file,timer,"Construction of Start System");
    hocotime := Elapsed_User_Time(timer);
  end Construct_Start_System;

  procedure Main is

    d,bz,bs,mv : natural;
    z : partition(p'range);
    nz : natural;
    mix : Integer_Vectors.Link_to_Vector;
    mixsub : Mixed_Subdivision;
    lifsup : Link_to_Array_of_Lists;
    qq : Poly_Sys(p'range);
    qqsols : Solution_List;

  begin
    Count_Roots(file,p,d,bz,bs,mv,z,nz,lifsup,mix,mixsub);
    Construct_Start_System
      (file,p,d,bz,bs,mv,z(1..nz),mix.all,lifsup.all,mixsub,rc,qq,qqsols);
    q := qq; qsols := qqsols;
    Clear(z);
    Clear(mix);
    Deep_Clear(lifsup);
    Deep_Clear(mixsub);
  end Main;

begin
  Main;
end Black_Box_Root_Counting;
