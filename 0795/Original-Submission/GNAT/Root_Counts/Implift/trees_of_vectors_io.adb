with integer_io;              use integer_io;
with Integer_Vectors;         use Integer_Vectors;
with Integer_Vectors_io;      use Integer_Vectors_io;

package body Trees_of_Vectors_io is

-- INTERNAL STATES :

  max : constant natural := 20;
  tokens : Vector(1..max);
  done : boolean;
  cnt : natural;

-- READ OPERATIONS :

  procedure get ( n : in natural; tv : in out Tree_of_Vectors ) is
  begin
    if done
     then cnt := 0;
          while not END_OF_LINE loop
            cnt := cnt + 1;
            get(tokens(cnt));
          end loop;
	  if cnt = 0
	   then done := true;
	   else done := false;
          end if;
    end if;
    if not done
     then if (cnt = 0) or (cnt > n)
	   then null;
	   elsif cnt = n
  	       then declare
	              nd : node;
                    begin
	              nd.d := new Integer_Vectors.Vector'(tokens(1..cnt));
		      done := true;
		      skip_line;
	              nd.ltv := new Tree_of_Vectors;
	              get(n-1,nd.ltv.all);
	              if Is_Null(nd.ltv.all)
	               then Clear(nd.ltv);
	                    nd.ltv := null;
                      end if;
	              Construct(nd,tv);
	              get(n,tv);
                    end;
               else get(n-1,tv);
          end if;
    end if;
  end get;

  procedure get ( file : in file_type;
		   n : in natural; tv : in out Tree_of_Vectors ) is
  begin
    if done
     then cnt := 0;
          while not END_OF_LINE(file) loop
            cnt := cnt + 1;
            get(file,tokens(cnt));
          end loop;
	  if cnt = 0
	   then done := true;
	   else done := false;
          end if;
    end if;
    if not done
     then if (cnt = 0) or (cnt > n)
	   then null;
	   elsif cnt = n
  	       then declare
	              nd : node;
                    begin
	              nd.d := new Integer_Vectors.Vector'(tokens(1..cnt));
		      done := true;
		      skip_line(file);
	              nd.ltv := new Tree_of_Vectors;
	              get(file,n-1,nd.ltv.all);
	              if Is_Null(nd.ltv.all)
	               then Clear(nd.ltv);
	                    nd.ltv := null;
                      end if;
	              Construct(nd,tv);
	              get(file,n,tv);
                    end;
               else get(file,n-1,tv);
          end if;
    end if;
  end get;

-- WRITE OPERATIONS :

  procedure put ( tv : in Tree_of_Vectors ) is
  begin
    put(Standard_Output,tv);
    new_line;
  end put;

  procedure put2 ( file : in file_type; tv : in Tree_of_Vectors );

  procedure put ( file : in file_type; nd : in node ) is
  begin
    put(file,nd.d); new_line(file);
    if not (nd.ltv = null) and then not Is_Null(nd.ltv.all)
     then put2(file,nd.ltv.all);
    end if;
  end put;

  procedure put2 ( file : in file_type; tv : in Tree_of_Vectors ) is

    tmp : Tree_of_Vectors := tv;

    procedure put_Node ( nd : in node; cont : out boolean ) is
    begin
      put(file,nd.d); new_line(file);
      if not (nd.ltv = null) and then not Is_Null(nd.ltv.all)
       then put2(file,nd.ltv.all);
      end if;
      cont := true;
    end put_Node;

   -- procedure put_Nodes is new Iterator ( Process => put_Node );

  begin
   -- put_Nodes(tv);
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put2;

  procedure put ( file : in file_type; tv : in Tree_of_Vectors ) is
  begin
    put2(file,tv);
    new_line(file);
  end put;

-- PACKAGE INITIALIZATION :

begin
  done := true;
end Trees_of_Vectors_io;
