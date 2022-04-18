with Face_Enumerators;                 use Face_Enumerators;

package body Integer_Faces_of_Polytope is

-- AUXILIAIRIES :

  function Create_Edge ( pts : Integer_Vectors_of_Vectors.Vector;
                         i,j : integer ) return Face is

  -- DESCRIPTION :
  --   Creates the edge spanned by pts(i) and pts(j).

    res : Face(0..1) := new Integer_Vectors_of_Vectors.Vector(0..1);

  begin
    res(0) := new Vector'(pts(i).all);
    res(1) := new Vector'(pts(j).all);
    return res;
  end Create_Edge;

  function Create_Face ( pts : Integer_Vectors_of_Vectors.Vector;
                         f : Vector ) return Face is

  -- DESCRIPTION :
  --   Returns vector of points pts(f(i)) that span the face.

    res : Face(f'range) := new Integer_Vectors_of_Vectors.Vector(f'range);

  begin
    for i in f'range loop
      res(i) := new Vector'(pts(f(i)).all);
    end loop;
    return res;
  end Create_Face;

  procedure Move_to_Front ( pts : in out Integer_Vectors_of_Vectors.Vector;
                            x : in Vector ) is

  -- DESCRIPTION :
  --   The vector x is move to the front of the vector pts.

  begin
    if pts(pts'first).all /= x
     then for i in pts'first+1..pts'last loop
            if pts(i).all = x
             then pts(i).all := pts(pts'first).all;
                  pts(pts'first).all := x;
                  return;
            end if;
          end loop;
    end if;
  end Move_to_Front;

-- CONSTRUCTORS :

  function Create ( k,n : positive; p : List ) return Faces is

    res : Faces;

  begin
    if k > n
     then return res;
     else
       declare
         m : constant natural := Length_Of(p);
         pts : Integer_Vectors_of_Vectors.Vector(1..m) := Shallow_Create(p);
         res_last : Faces := res;
       begin
         if k = 1
          then
            declare
              procedure Append_Edge ( i,j : in natural; cont : out boolean ) is
                f : Face := Create_Edge(pts,i,j);
              begin
                Append(res,res_last,f); cont := true;
              end Append_Edge;
              procedure Enum_Edges is new Enumerate_Edges(Append_Edge);
            begin
              Enum_Edges(pts);
            end;
          else
            declare
              procedure Append_Face ( fa : in Vector; cont : out boolean ) is
                f : Face := Create_Face(pts,fa);
              begin
                Append(res,res_last,f); cont := true;
              end Append_Face;
              procedure Enum_Faces is new Enumerate_Faces(Append_Face); 
            begin
              Enum_Faces(k,pts);
            end;
         end if;
         return res;
       end;
    end if;
  end Create;

  function Create ( k,n : positive; p : List; x : Vector ) return Faces is

    res : Faces;

  begin
    if k > n
     then return res;
     else
       declare
         m : constant natural := Length_Of(p);
         pts : Integer_Vectors_of_Vectors.Vector(1..m) := Shallow_Create(p);
         res_last : Faces := res;
       begin
         Move_to_Front(pts,x);
         if k = 1
          then
            declare
              procedure Append_Edge ( i,j : in natural; cont : out boolean ) is
                f : Face;
              begin
                if i = pts'first
                 then f := Create_Edge(pts,i,j);
                      Append(res,res_last,f);
                      cont := true;
                 else cont := false;
                end if;
              end Append_Edge;
              procedure Enum_Edges is new Enumerate_Edges(Append_Edge);
            begin
              Enum_Edges(pts);
            end;
          else
            declare
              procedure Append_Face ( fa : in Vector; cont : out boolean ) is
                f : Face;
              begin
                if fa(fa'first) = pts'first
                 then f := Create_Face(pts,fa);
                      Append(res,res_last,f);
                      cont := true;
                 else cont := false;
                end if;
              end Append_Face;
              procedure Enum_Faces is new Enumerate_Faces(Append_Face);
            begin
              Enum_Faces(k,pts);
            end;
         end if;
         return res;
       end;
    end if;
  end Create;

  function Create_Lower ( k,n : positive; p : List ) return Faces is

    res : Faces;

  begin
    if k > n
     then return res;
     else
       declare
         m : constant natural := Length_Of(p);
         pts : Integer_Vectors_of_Vectors.Vector(1..m) := Shallow_Create(p);
         res_last : Faces := res;
       begin
         if k = 1
          then
            declare
              procedure Append_Edge ( i,j : in natural; cont : out boolean ) is
                f : Face := Create_Edge(pts,i,j);
              begin
                Append(res,res_last,f); cont := true;
              end Append_Edge;
              procedure Enum_Edges is new Enumerate_Lower_Edges(Append_Edge);
            begin
              Enum_Edges(pts);
            end;
          else
            declare
              procedure Append_Face ( fa : in Vector; cont : out boolean ) is
                f : Face := Create_Face(pts,fa);
              begin
                Append(res,res_last,f); cont := true;
              end Append_Face;
              procedure Enum_Faces is new Enumerate_Lower_Faces(Append_Face);
            begin
              Enum_Faces(k,pts);
            end;
         end if;
         return res;
       end;
    end if;
  end Create_Lower;

  function Create_Lower ( k,n : positive; p : List; x : Vector )
                        return Faces is

    res : Faces;

  begin
    if k > n
     then return res;
     else
       declare
         m : constant natural := Length_Of(p);
         pts : Integer_Vectors_of_Vectors.Vector(1..m) := Shallow_Create(p);
         res_last : Faces := res;
       begin
         Move_to_Front(pts,x);
         if k = 1
          then
            declare
              procedure Append_Edge ( i,j : in natural; cont : out boolean ) is
                f : Face := Create_Edge(pts,i,j);
              begin
                if i = pts'first
                 then f := Create_Edge(pts,i,j);
                      Append(res,res_last,f);
                      cont := true;
                 else cont := false;
                end if;
              end Append_Edge;
              procedure Enum_Edges is new Enumerate_Lower_Edges(Append_Edge);
            begin
              Enum_Edges(pts);
            end;
          else
            declare
              procedure Append_Face ( fa : in Vector; cont : out boolean ) is
                f : Face;
              begin
                if fa(fa'first) = pts'first
                 then f := Create_Face(pts,fa);
                      Append(res,res_last,f);
                      cont := true;
                 else cont := false;
                end if;
              end Append_Face;
              procedure Enum_Faces is new Enumerate_Lower_Faces(Append_Face);
            begin
              Enum_Faces(k,pts);
            end;
         end if;
         return res;
       end;
    end if;
  end Create_Lower;

  procedure Construct ( first : in out Faces; fs : in Faces ) is

    tmp : Faces := fs;

  begin
    while not Is_Null(tmp) loop
      Construct(Head_Of(tmp),first);
      tmp := Tail_Of(tmp);
    end loop;
  end Construct;

  procedure Copy ( f1 : in Face; f2 : in out Face ) is
  begin
    Deep_Clear(f2);
    f2 := new Integer_Vectors_of_Vectors.Vector(f2'range);
    for i in f2'range loop
      f2(i) := new Integer_Vectors.Vector'(f1(i).all);
    end loop;
  end Copy;

  procedure Deep_Copy ( f1 : in Faces; f2 : in out Faces ) is

    tmp,last : Faces;

  begin
    Deep_Clear(f2);
    tmp := f1;
    while not Is_Null(tmp) loop
      declare
        face1 : Face := Head_Of(tmp);
        face2 : Face := new Integer_Vectors_of_Vectors.Vector(face1'range);
      begin
        for i in face2'range loop
          face2(i) := new Integer_Vectors.Vector'(face1(i).all);
        end loop;
        Append(f2,last,face2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Deep_Copy;

-- SELECTORS :

  function Is_Equal ( f1,f2 : Face ) return boolean is

    found : boolean;

  begin
    for i in f1'range loop
      found := false;
      for j in f2'range loop
        found := Equal(f1(i).all,f2(j).all);
        exit when found;
      end loop;
      if not found
       then return false;
      end if;
    end loop;
    return true;
  end Is_Equal;

  function Is_In ( f : Face; x : Vector ) return boolean is
  begin
    for i in f'range loop
      if f(i).all = x
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( fs : Faces; f : Face ) return boolean is

    tmp : Faces := fs;

  begin
    while not Is_Null(tmp) loop
      if Is_Equal(f,Head_Of(tmp))
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Extract_Faces ( fs : Faces; x : Vector ) return Faces is

    res,res_last,tmp : Faces;

  begin
    tmp := fs;
    while not Is_Null(tmp) loop
      declare
        f : Face := Head_Of(tmp);
      begin
        if Is_In(f,x)
         then Append(res,res_last,f);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Extract_Faces;

-- DESTRUCTORS :

  procedure Deep_Clear ( f : in out Face ) is
  begin
    if f /= null
     then for i in f'range loop
            Clear(f(i));
          end loop;
    end if;
    Shallow_Clear(f);
  end Deep_Clear;

  procedure Shallow_Clear ( f : in out Face ) is
  begin
    if f /= null
     then Integer_Vectors_of_Vectors.Clear(f.all);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( fa : in out Face_Array ) is
  begin
    for i in fa'range loop
      Deep_Clear(fa(i));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( fa : in out Face_Array ) is
  begin
    for i in fa'range loop
      Shallow_Clear(fa(i));
    end loop;
  end Shallow_Clear;

  procedure Deep_Clear ( fs : in out Faces ) is

    tmp : Faces := fs;

  begin
    while not Is_Null(tmp) loop
      declare
	f : Face := Head_Of(tmp);
      begin
	Deep_Clear(f);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Faces.Clear(Lists_of_Faces.List(fs));
  end Deep_Clear;

  procedure Shallow_Clear ( fs : in out Faces ) is

    tmp : Faces := fs;

  begin
    Lists_of_Faces.Clear(Lists_of_Faces.List(fs));
  end Shallow_Clear;

  procedure Deep_Clear ( afs : in out Array_of_Faces ) is
  begin
    for i in afs'range loop
      Deep_Clear(afs(i));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( afs : in out Array_of_Faces ) is
  begin
    for i in afs'range loop
      Shallow_Clear(afs(i));
    end loop;
  end Shallow_Clear;

end Integer_Faces_of_Polytope;
